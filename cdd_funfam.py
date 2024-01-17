#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script generates a list of CDD entries overlapping with CATH-Funfams

@arguments [CONFIG_FILE]: file containing database connection and files location (see config.ini)

"""

import argparse
import os
import re
from configparser import ConfigParser
import traceback
import cx_Oracle
from statistics import mean 
import go_comparison.funfam2gomapping as f2go

def get_overlap(start1, end1, start2, end2):
    totalRange = max(end1, end2) - min(start1, start2)
    sumOfRanges = (end1 - start1) + (end2 - start2)
    overlappingInterval = 0

    if sumOfRanges > totalRange: # means they overlap
        overlappingInterval = min(end1, end2) - max(start1, start2)

    return overlappingInterval

def get_total_cdd_matches(connectString):
    try:
        # subprocess.run(["source", "~oracle/ora112setup.sh"])
        con = cx_Oracle.connect(connectString)
        cur = con.cursor()
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)

    request = """
        select method_ac, count(distinct protein_ac)
        from match
        where dbcode='J'
        group by method_ac
    """
    cur.execute(request)
    cdd_matches = dict(cur.fetchall())

    cur.close()
    con.close()

    return cdd_matches


def find_overlaping(user, password, schema, outputdir, outfile_one_match, outfile_multi_match, funfam2gofile):
    connectString = "".join([user, "/", password, "@", schema])
    list_cdd = dict()
    list_funfam = set()

    cdd_matches = get_total_cdd_matches(connectString)

    try:
        # subprocess.run(["source", "~oracle/ora112setup.sh"])
        connection = cx_Oracle.connect(connectString)
        cursor = connection.cursor()
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)

    #get funfam2go list from file
    funfam2go = dict()
    if os.path.exists(funfam2gofile) and os.path.getsize(funfam2gofile) > 0:
        funfam2go = f2go.get_mapping_from_file(funfam2gofile)

    request="""
        SELECT m.protein_ac, m.method_ac, m.pos_from, m.pos_to, fm.method_ac, fm.pos_from, fm.pos_to
        FROM INTERPRO.MATCH m
        JOIN INTERPRO.FEATURE_MATCH fm on m.protein_ac=fm.protein_ac
        LEFT JOIN INTERPRO.ENTRY2METHOD e2m on m.method_ac=e2m.method_ac
        LEFT JOIN INTERPRO.METHOD_COMMENT m2c on m2c.method_ac=e2m.method_ac
        WHERE m.DBCODE='J' AND fm.dbcode='f'
        AND e2m.method_ac is null
        AND m2c.method_ac is null
    """

    cursor.execute(request)
    for row in cursor:
        protein, cdd, cdd_start, cdd_end, funfam, funfam_start, funfam_end = row

        cdd_len = int(cdd_end)-int(cdd_start)
        funfam_len = int(funfam_end)-int(funfam_start)
        len_overlap = get_overlap(int(cdd_start), int(cdd_end), int(funfam_start), int(funfam_end))

        if funfam not in funfam2go:
            list_funfam.add((outputdir, funfam))

        try:
            list_cdd[cdd][funfam]["protein"].append(protein)
            list_cdd[cdd][funfam]["cdd_doms"].append(cdd_len)
            list_cdd[cdd][funfam]["funfam_doms"].append(funfam_len)
        except KeyError:
            try:
                list_cdd[cdd][funfam] = {"protein":[protein], "overlap":[], "cdd_doms":[cdd_len], "funfam_doms":[funfam_len]}
            except KeyError:
                list_cdd[cdd] = {funfam: {"protein":[protein], "overlap":[], "cdd_doms":[cdd_len], "funfam_doms":[funfam_len]}}
        
        if len_overlap > cdd_len/2 and len_overlap > funfam_len/2: #overlap should be at least half the length of the CDD and funfam domains
            list_cdd[cdd][funfam]["overlap"].append(len_overlap)
            
    cursor.close()
    connection.close()

    #search funfam GO terms and save results in file for subsequent runs
    if len(list_funfam) > 0:
        funfam2go.update(f2go.get_list_funfam2go(list(list_funfam)))
        f2go.save_in_file(funfam2gofile, funfam2go)

    with open(outfile_one_match, "w") as f1, open(outfile_multi_match, "w") as fm:
        header = "cdd; link to cdd website; funfam; total number of CDD/UniProt matches; total number of UniProt matches for the CDD/funfam pair; number of UniProt matches that overlap for the CDD/funfam pair;  percentage of CDD/UniProt matches vs overlap; percentage of proteins overlapping; average overlap length of CDD/funfam; average length of the CDD domain; average length of the funfam domain; comments\n"
        # header = "cdd, funfam, total number of CDD/UniProt matches, total number of UniProt matches for the CDD/funfam pair, number of UniProt matches that overlap for the CDD/funfam pair, percentage of CDD/UniProt matches vs overlap, percentage of proteins overlapping, average overlap length of CDD/funfam, average length of the CDD domain, average length of the funfam domain, comments\n"

        f1.write(header)
        fm.write(header)

        for cdd, funfams in list_cdd.items():
            count_overlap = 0
            text = ""
            
            prontolink = f'=HYPERLINK("http://pronto.ebi.ac.uk:5000/signature/{cdd}","{cdd}")'
            cddlink = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={cdd}","ncbi_{cdd}")'
            for funfam, content in funfams.items():
                percent_cdd_match_overlap = round(len(content['overlap'])*100/int(cdd_matches[cdd]), 2)
                percentoverlap = round(len(content['overlap'])*100/len(content['protein']), 2)
                nb_funfam2go = len(funfam2go[funfam]) if funfam in funfam2go else 0
                if len(content["overlap"]) > 0 and percentoverlap >= 50 and percent_cdd_match_overlap >= 50 and nb_funfam2go > 0 and nb_funfam2go < 5:
                    count_overlap += 1 
                    gene3d = re.match(r'G3DSA:((\d+\.)+\d+):FF:0+([1-9]\d*)', funfam).group(1)
                    funfamid = re.match(r'G3DSA:((\d+\.)+\d+):FF:0+([1-9]\d*)', funfam).group(3)
                    funfamlink = f'=HYPERLINK("http://www.cathdb.info/version/latest/superfamily/{gene3d}/funfam/{funfamid}","{funfam}")'
                    
                    text+=f"{prontolink};{cddlink};{funfamlink};{cdd_matches[cdd]};{len(content['protein'])};{len(content['overlap'])};{percent_cdd_match_overlap};{percentoverlap};{round(mean(content['overlap']),2)};{round(mean(content['cdd_doms']),2)};{round(mean(content['funfam_doms']),2)};\n"
                    # text+=f"{cdd},{funfam},{percent_cdd_match_overlap},{len(content['protein'])},{len(content['overlap'])},{percent_cdd_match_overlap},{percentoverlap},{round(mean(content['overlap']),2)},{round(mean(content['cdd_doms']),2)},{round(mean(content['funfam_doms']),2)},\n"

            if count_overlap == 1:
                f1.write(text)
            else:
                fm.write(text)
                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    outputdir = config["dir"]["outputdir"]
    outfile_one_match = os.path.join(outputdir, "list_cdd_funfam_one_match_test.csv")
    outfile_multi_match = os.path.join(outputdir, "list_cdd_funfam_multi_match_test.csv")
    funfam2gofile = os.path.join(outputdir, "funfam2go.csv")

    # database connection values
    user = config["database"]["user"]
    password = config["database"]["password"]
    schema = config["database"]["schema"]

    find_overlaping(user, password, schema, outputdir, outfile_one_match, outfile_multi_match, funfam2gofile)