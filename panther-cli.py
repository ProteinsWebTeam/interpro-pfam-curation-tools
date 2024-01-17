"""
# Find potential replacements for deleted signatures when updating PANTHER

# Export PANTHER matches (memory: ~2GB)
$ export INTERPRO_URL="interpro/*******@IPPRO"
$ python panther-cli.py export -a 84 -o /hps/nobackup/agb/interpro/typhaine/panther/panther15
$ python panther-cli.py export -a 101 -o /hps/nobackup/agb/interpro/typhaine/panther/panther17

# Find replacements for deleted signatures (memory: ~20GB)
$ python panther-cli.py find -1 /hps/nobackup/agb/interpro/mblum/panther/panther15 -2 /hps/nobackup/agb/interpro/typhaine/panther/panther17 > /hps/nobackup/agb/interpro/typhaine/panther/replacements.tsv

# List deleted signatures
$ python panther-cli.py list -1 /hps/nobackup/agb/interpro/mblum/panther/panther15 -2 /hps/nobackup/agb/interpro/typhaine/panther/panther17

# List gained and lost proteins (memory: ~8GB)
$ python panther-cli.py diff -1 /hps/nobackup/agb/interpro/mblum/panther/panther15 -2 /hps/nobackup/agb/interpro/typhaine/panther/panther17 > /hps/nobackup/agb/interpro/typhaine/panther/sequences.tsv
"""

import argparse
import heapq
import os
import pickle
import re
import shutil
import struct
import sys
import tempfile
import psycopg2

import cx_Oracle


MIN_SIMILARITY = 0.80


class File:
    def __init__(self, path, cached=False):
        self.path = path
        self.fh = None
        self.footer = {}
        self.cached = cached
        self.data = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __contains__(self, item):
        return item in self.footer

    def __iter__(self):
        for key in sorted(self.footer):
            yield key

    def __getitem__(self, key):
        try:
            return self.data[key]
        except KeyError:
            pass

        offset = self.footer[key]
        self.fh.seek(offset)
        val = ProteinSet(pickle.load(self.fh))
        if self.cached:
            self.data[key] = val
        return val

    def __len__(self):
        return len(self.footer)

    @property
    def signatures(self):
        return list(self.footer.keys())

    def write(self, files):
        with open(self.path, "wb") as fh:
            offset = fh.write(struct.pack("<Q", 0))
            footer = {}

            iterables = [self.load(path) for path in files]
            acc = None
            proteins = []

            for key, values in heapq.merge(*iterables):
                if key != acc:
                    if acc:
                        footer[acc] = offset
                        offset += fh.write(pickle.dumps(proteins))

                    acc = key
                    proteins.clear()

                proteins += values

            if acc:
                footer[acc] = offset
                offset += fh.write(pickle.dumps(proteins))

            pickle.dump(footer, fh)
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def open(self):
        self.close()
        self.fh = open(self.path, "rb")
        offset, = struct.unpack("<Q", self.fh.read(8))
        self.fh.seek(offset)
        self.footer = pickle.load(self.fh)

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

    @staticmethod
    def dump(cache, file):
        with open(file, "wb") as fh:
            for key in sorted(cache):
                pickle.dump((key, cache[key]), fh)

    @staticmethod
    def load(file):
        with open(file, "rb") as fh:
            while True:
                try:
                    item = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield item


class ProteinSet:
    def __init__(self, proteins):
        self.data = proteins

    @property
    def all(self):
        return {acc for acc, reviewed in self.data}

    @property
    def reviewed(self):
        return {acc for acc, reviewed in self.data if reviewed}

def get_fragments(pronto_uri):
    con = connect_pg(pronto_uri)
    cur = con.cursor()

    cur.execute("""
        SELECT accession 
        FROM protein
        WHERE is_fragment = 't'
        AND is_reviewed = 't'
    """)

    fragments = set(cur.fetchall())
    cur.close()
    con.close()

    return fragments


def export(uri, pronto_uri, args):
    
    file = File(args.o)
    tmpdir = f"{args.o}_tmp"
    try:
        shutil.rmtree(tmpdir)
    except FileNotFoundError:
        pass

    os.makedirs(tmpdir)

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT M.METHOD_AC, P.AC, P.DBID
        FROM IPRSCAN.IPM_PANTHER_MATCH@ISPRO M
        INNER JOIN (
            SELECT UPI, AC, DBID
            FROM UNIPARC.XREF
            WHERE DBID = 2
            AND DELETED = 'N'
        ) P ON M.UPI = P.UPI
        WHERE M.ANALYSIS_ID = :1
        """,
        [args.a]
    )

    cache = {}
    files = []
    i = 0
    for method_acc, protein_acc, database_id in cur:
        is_reviewed = database_id == 2

        # if protein_acc in fragments:
        #     continue
        
        try:
            cache[method_acc].append((protein_acc, is_reviewed))
        except KeyError:
            cache[method_acc] = [(protein_acc, is_reviewed)]

        i += 1
        if i % 1e6 == 0:
            fd, path = tempfile.mkstemp(dir=tmpdir)
            file.dump(cache, fd)
            cache.clear()
            files.append(path)

    cur.close()
    con.close()

    if cache:
        fd, path = tempfile.mkstemp(dir=tmpdir)
        file.dump(cache, fd)
        cache.clear()
        files.append(path)

    file.write(files)
    shutil.rmtree(tmpdir)

def connect_pg(url):
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', url)
    return psycopg2.connect(
        user=m.group(1),
        password=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        dbname=m.group(5)
    )

def get_integrated(url):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M ON EM.METHOD_AC = M.METHOD_AC
        WHERE M.DBCODE = 'V'
        """
    )
    integrated = dict(cur.fetchall())
    cur.close()
    con.close()
    return integrated

def get_deleted_from_file(del_file):
    deleted = dict()
    with open (del_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            entry_ac, method_ac = line.split('\t')
            deleted[method_ac] = entry_ac

    return deleted

def get_integrated_with_otherdb(url):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M ON EM.METHOD_AC = M.METHOD_AC
        WHERE M.DBCODE != 'V'
        """
    )
    integrated = set(cur.fetchall())
    cur.close()
    con.close()
    return integrated

def list_deleted(url, pronto_url, args):
    integrated = get_integrated(url)

    with File(args.f1) as now, File(args.f2) as nxt:
        now.open()
        nxt.open()
        for s_acc in sorted(integrated.keys()):
            if s_acc in now and s_acc not in nxt:
                print(f"{s_acc}\t{integrated[s_acc]}")


def gen_flag(similarity, same_reviewed):
    if similarity >= 0.9 and same_reviewed:
        return f"substitute ({similarity*100:.0f}%)"
    elif same_reviewed:
        return f"candidate ({similarity*100:.0f}%)"
    else:
        return f"{similarity*100:.0f}%"


def replacements_otherdb(url, pronto_url):

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M ON EM.METHOD_AC = M.METHOD_AC
        WHERE M.DBCODE != 'V'
        """
    )
    signatures = dict(cur.fetchall())
    cur.close()
    con.close()

    con = connect_pg(pronto_url)
    cur = con.cursor()
    # print("search otherdb")

    cur.execute("""
        SELECT signature_acc, protein_acc
        FROM signature2protein
        WHERE is_reviewed = 't' 
        AND signature_acc NOT LIKE 'PTHR%'
    """)

    replacement = dict()
    for row in cur.fetchall():
        sign = row[0]
        prot = row[1]
        try: 
            replacement[sign].add(prot)
        except KeyError:
            replacement[sign] = set([prot])
    
    cur.close()
    con.close()

    # print(len(replacement))
    return replacement, signatures


def find_replacements(uri, pronto_uri, args):
    fragments = get_fragments(pronto_uri)
    # integrated = get_integrated(uri)
    integrated = get_deleted_from_file(args.f3)

    # entries_with_other_db = get_integrated_with_otherdb(uri)
    other_db, signatures = replacements_otherdb(uri, pronto_uri)
    count = 0 

    with File(args.f1) as now, File(args.f2, cached=True) as nxt:
        now.open()
        nxt.open()
        for s_acc in sorted(integrated):
            if s_acc not in now or s_acc in nxt:
                continue            
            
            e_acc = integrated[s_acc]
            proteins = now[s_acc]
            
            reviewed = proteins.reviewed.difference(fragments)
            # print(s_acc, e_acc, reviewed)
            
            # if e_acc in entries_with_other_db:
            #     print(f"{s_acc}\t{e_acc}\tintegrated with other signatures\t-\t-")
            #     continue

            m = re.fullmatch(r"(PTHR\d+):SF\d+", s_acc)
            parent_acc = m.group(1) if m else None

            candidates = []
            for other_acc in nxt:
                other_proteins = nxt[other_acc]
                op_reviewed = other_proteins.all.difference(fragments)
                common = len(reviewed & op_reviewed)
                if common == 0:
                    continue

                union = len(reviewed) + len(op_reviewed) - common
                similarity = common / union

                if other_acc == parent_acc or similarity >= MIN_SIMILARITY:
                    candidates.append((
                        other_acc,
                        integrated.get(other_acc, ''),
                        similarity,
                        reviewed and reviewed == op_reviewed
                    ))

            for sign_acc in other_db:
                common = len(reviewed & other_db[sign_acc])
                if common == 0:
                    continue

                union = len(reviewed) + len(other_db[sign_acc]) - common
                similarity = common / union

                if similarity >= MIN_SIMILARITY:
                    candidates.append((
                        sign_acc,
                        signatures.get(sign_acc, ''),
                        similarity,
                        reviewed and reviewed == other_db[sign_acc]
                    ))

            if not candidates:
                print(f"{s_acc}\t{e_acc}\t-\t-\t-")
                continue

            candidates.sort(key=lambda x: -x[2])

            other_acc, other_entry, similarity, reviewed = candidates[0]
            flag = gen_flag(similarity, reviewed)
            print(f"{s_acc}\t{e_acc}\t{other_acc}\t{other_entry}\t{flag}")
            count+=1

            for other_acc, other_entry, similarity, reviewed in candidates[1:]:
                flag = gen_flag(similarity, reviewed)
                print(f"\t\t{other_acc}\t{other_entry}\t{flag}")

            # if count >= 100:
            #     break

def find_sequences(uri, pronto_uri, args):
    integrated = get_integrated(uri)

    rev_now = set()
    unr_now = set()
    with File(args.f1) as f:
        f.open()
        for s_acc in sorted(integrated.keys()):
            try:
                proteins = f[s_acc]
            except KeyError:
                continue
            else:
                for acc, reviewed in proteins.data:
                    if reviewed:
                        rev_now.add(acc)
                    else:
                        unr_now.add(acc)

    rev_gained = []
    unr_gained = []
    with File(args.f2) as f:
        f.open()
        for s_acc in sorted(integrated.keys()):
            try:
                proteins = f[s_acc]
            except KeyError:
                continue
            else:
                for acc, reviewed in proteins.data:
                    if reviewed:
                        if acc in rev_now:
                            rev_now.remove(acc)
                        else:
                            rev_gained.append(acc)
                    elif acc in unr_now:
                        unr_now.remove(acc)
                    else:
                        unr_gained.append(acc)

    for acc in sorted(rev_gained):
        print(f"{acc}\tswissprot\tgained")

    for acc in sorted(unr_gained):
        print(f"{acc}\ttrembl\tgained")

    for acc in sorted(rev_now):
        print(f"{acc}\tswissprot\tlost")

    for acc in sorted(unr_now):
        print(f"{acc}\ttrembl\tlost")


def main():
    parser = argparse.ArgumentParser(description="PANTHER update helper")
    subparsers = parser.add_subparsers()

    parser_exp = subparsers.add_parser("export", help="export files")
    parser_exp.add_argument("-a", help="analysis ID", type=int, required=True)
    parser_exp.add_argument("-o", help="output file", required=True)
    parser_exp.set_defaults(func=export)

    parser_del = subparsers.add_parser("list",
                                       help="list deleted signatures")
    parser_del.add_argument("-1", dest="f1", help="current version file",
                            required=True)
    parser_del.add_argument("-2", dest="f2", help="next version file",
                            required=True)
    parser_del.set_defaults(func=list_deleted)

    parser_rep = subparsers.add_parser("find",
                                        help="find replacements "
                                             "for deleted signatures")
    parser_rep.add_argument("-1", dest="f1", help="current version file",
                            required=True)
    parser_rep.add_argument("-2", dest="f2", help="next version file",
                            required=True)
    parser_rep.add_argument("-3", dest="f3", help="deleted signatures file",
                            required=True)
    parser_rep.set_defaults(func=find_replacements)

    parser_rep = subparsers.add_parser("diff",
                                        help="find gained/lost sequences")
    parser_rep.add_argument("-1", dest="f1", help="current version file",
                            required=True)
    parser_rep.add_argument("-2", dest="f2", help="next version file",
                            required=True)
    parser_rep.set_defaults(func=find_sequences)

    args = parser.parse_args()

    try:
        uri = os.environ["INTERPRO_URL"]
    except KeyError:
        parser.error("Environment variable INTERPRO_URL not defined")
    try:
        pronto_uri = os.environ["PRONTO_URL"]
    except KeyError:
        parser.error("Environment variable PRONTO_URL not defined")
    print(args.func)
#    args.func(uri, pronto_uri, args)


if __name__ == '__main__':
    main()
