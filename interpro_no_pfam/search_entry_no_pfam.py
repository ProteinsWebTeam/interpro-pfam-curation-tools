import argparse
import sys, os, re
from configparser import ConfigParser
import cx_Oracle
import psycopg
import traceback
from multiprocessing import Pool

def getOracleConnection(connectString):
    """
    Set database connection
    """
    try:
        connection = cx_Oracle.connect(connectString)
        cursor = connection.cursor()

        return connection, cursor
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)
        if "invalid username" in stackTrace:
            print("Could not connect to {0} as user {1}".format(schema, user))
            print(
                "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
            )
            print('eg "ops$craigm" or ops\$craigm')
            print("Otherwise the shell will remove the '$' and all subsequent characters!")
            sys.exit(1)


def getPostgresConnection(connectString):
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', connectString)

    if m is None:
        raise RuntimeError(f"invalid connection string: {connectString}")

    return psycopg.connect(
        user=m.group(1),
        password=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        dbname=m.group(5)
    )


def split_list(l, n):
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i : i + n]


def get_proteins_for_signatures(cursor, signatures):
    placeholders = ', '.join(['%s'] * len(signatures))
    cursor.execute(f"SELECT PROTEIN_ACC FROM SIGNATURE2PROTEIN WHERE SIGNATURE_ACC IN ({placeholders})", (*signatures,))
    protein_list = set(row[0] for row in cursor)

    return protein_list


def get_proteins_no_pfam(cursor, proteins):
    all_proteins_pfam = set()
    if len(proteins) > 1000:
        chunks = split_list(list(proteins), 1000)
        for chunk in chunks:
            placeholders = ', '.join(['%s'] * len(chunk))
            cursor.execute(f"SELECT PROTEIN_ACC FROM SIGNATURE2PROTEIN WHERE PROTEIN_ACC IN ({placeholders}) AND SIGNATURE_ACC LIKE %s", (*chunk, 'PF%'))
            protein_list = set(row[0] for row in cursor)
            all_proteins_pfam = all_proteins_pfam.union(protein_list)   
    else:
        placeholders = ', '.join(['%s'] * len(proteins))
        cursor.execute(f"SELECT PROTEIN_ACC FROM SIGNATURE2PROTEIN WHERE PROTEIN_ACC IN ({placeholders}) AND SIGNATURE_ACC LIKE %s", (*proteins, 'PF%'))
        protein_list = set(row[0] for row in cursor)
        all_proteins_pfam = all_proteins_pfam.union(protein_list)
    
    return all_proteins_pfam


def get_proteins_for_pfam(cursor):
    print("Getting proteins for Pfam")
    cursor.execute(f"SELECT PROTEIN_ACC FROM SIGNATURE2PROTEIN WHERE SIGNATURE_ACC LIKE 'PF%'")
    protein_list = set(row[0] for row in cursor)
    print(len(protein_list))
    return protein_list


def exclude_entries_with_pfam(cursor, entry_ac):
    print("Excluding entries with pfam signatures")
    cursor.execute(f"SELECT METHOD_AC FROM INTERPRO.ENTRY2METHOD WHERE ENTRY_AC=:1", (entry_ac,))

    entry2method = list()
    for item in cursor.fetchall():
        if re.search("^PF", item[0]):
            return []
        entry2method.append(item[0])
    return entry2method


def get_entries(oracle_connectionString, postgresql_connectionString, inputfile, processed_entries, outputdir):
    connection, cursor = getOracleConnection(oracle_connectionString)
    entries_excluded = set()
    entries_processed = set()
    entries = []

    with open(inputfile, "r") as f:
        for line in f.readlines():
            entry_ac = line.strip("\n")
            outputfile = os.path.join(outputdir, f"{entry_ac}.tsv")

            if entry_ac in processed_entries or os.path.exists(outputfile):
                entries_processed.add(entry_ac)
                continue

            signatures = exclude_entries_with_pfam(cursor, entry_ac)
            if len(exclude_entries_with_pfam(cursor, entry_ac)) == 0:
                entries_excluded.add(entry_ac)
                continue
            else:
                tuple = (entry_ac, postgresql_connectionString, signatures, outputfile)
                entries.append(tuple)

    cursor.close()
    connection.close()

    print(f"total number of entries: {len(entries_excluded) + len(entries) + len(entries_processed)}")
    print(f"number of entries excluded: {len(entries_excluded)}")
    print(f"number of entries already processed: {len(entries_processed)}")
    print(f"number of entries to process: {len(entries)}")

    return entries

def process_entry(entry_ac, postgresql_connectionString, signatures, outputfile):

    pg_con = getPostgresConnection(postgresql_connectionString)
    pg_cur = pg_con.cursor()

    print(entry_ac, signatures)
    sign_proteins = get_proteins_for_signatures(pg_cur, signatures)
    proteins_pfam = get_proteins_no_pfam(pg_cur, sign_proteins)
    proteins_no_pfam = sign_proteins - proteins_pfam

    with open(outputfile, "w") as outfile:
        outfile.write(f"{entry_ac}\t{len(sign_proteins)}\t{len(proteins_pfam)}\t{len(proteins_no_pfam)}\n")

    pg_cur.close()
    pg_con.close()

    return entry_ac

    
def get_proteins_for_entries(entries):

    with Pool(10) as p:
        results = p.starmap(process_entry, entries)

    print(len(results))


def find_processed_entries(processedfile):
    entries = []
    if os.path.exists(processedfile):
        with open(processedfile, "r") as f:
            for line in f.readlines():
                entry_ac = line.strip("\n").split("\t")[0]
                entries.append(entry_ac)

    return entries

            
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    oracle_connectionString = config["database"]["ipro-interpro"]
    postgresql_connectionString = config["database"]["pronto"]
    outputdir = config["files"]["outputdir"]
    inputfile = config["files"]["inputfile"]
    outputfile = config["files"]["outputfile"]
    processed_file = config["files"]["outputfile_processed"]
    
    # search if any entries have already been processed
    processed_entries = find_processed_entries(processed_file)

    # get the entries to search
    entries_to_search = get_entries(oracle_connectionString, postgresql_connectionString, inputfile, processed_entries, outputdir)
    
    #generate number of proteins and overlap with Pfam for each entry (outputdir/entry_id.tsv)
    get_proteins_for_entries(entries_to_search)

    # concatenate all the files in outputdir
    with open(outputfile, 'w') as outfile:
        os.chdir(outputdir)
        tsv_files = glob.glob("*.tsv")
        command = ["cat"] + tsv_files

        subprocess.run(command, stdout=outfile)
