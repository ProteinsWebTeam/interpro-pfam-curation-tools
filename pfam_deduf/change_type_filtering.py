import cx_Oracle
import traceback
import sys, errno, re, json, ssl, os
from urllib import request
from urllib.error import HTTPError
from time import sleep
import argparse
from configparser import ConfigParser
from multiprocessing import Pool

def db_connection(uri: str):
    try:
        connection = cx_Oracle.connect(uri)
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)

    return connection

def get_pfam_not_domain(uri):
    con = db_connection(uri)
    cur = con.cursor()

    request = """select m.method_ac,  m.sig_type, e.entry_ac, e.entry_type
                from InterPro.entry e
                join interpro.entry2method em on em.entry_ac=e.entry_ac
                join interpro.method m on m.method_ac=em.method_ac
                where e.entry_type='D' 
                and m.sig_type='F' 
                and m.method_ac like 'PF%'
    """
    # and m.method_ac like 'PF%'
    # and m.method_ac='PF00026'
    # and m.method_ac='PF03121'
    cur.execute(request)

    list_doms = {}

    for row in cur:
        list_doms[row[0]] = {"m_type":row[1], "e_acc":row[2], "e_type":row[3]}

    cur.close()
    con.close()

    return list_doms

def get_pfam_ida(accession: str):
    url=f"https://www.ebi.ac.uk/interpro/api/entry?ida_search={accession}&page_size=200"

    total_nb_prot = get_total_protein_count(accession)

    #disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next = url
  
    attempts = 0
    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                #no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next = payload["next"]
            attempts = 0
            if not next:
                last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
            # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

        for i, item in enumerate(payload["results"]):
            protein_count = item["unique_proteins"]
            domains = item["ida"].split('-')
            if len(domains) == 1 and protein_count > total_nb_prot/2:
                prot_length = item["representative"]["length"]
                protein = item["representative"]["accession"]
                start_dom = item["representative"]["domains"][0]["coordinates"][0]["fragments"][0]["start"]
                end_dom = item["representative"]["domains"][0]["coordinates"][0]["fragments"][0]["end"]
                dom_length = end_dom-start_dom
                if dom_length < prot_length/2: #pfam should cover at least 1/2 protein length
                    return domains, protein_count, total_nb_prot, protein, dom_length, prot_length
        
        # Don't overload the server, give it time before asking for more
        if next:
            sleep(1)

    return None

def get_total_protein_count(accession):
    url=f"https://www.ebi.ac.uk/interpro/api/protein/entry/pfam/{accession}"

    #disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next = url
  
    attempts = 0
    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                #no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next = ""
            attempts = 0
            if not next:
                last_page = True
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
            # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

        count = (payload["proteins"]["uniprot"])
    
    return count
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    # database connection uri
    uri = config["database"]["ipro-interpro"]

    domains_to_check = {}

    accessions = get_pfam_not_domain(uri)

    with Pool(10) as p:
        results = p.map(get_pfam_ida, accessions)

    for item in results:
        if item:
            ida, count_ida, count_total, protein, dom_length, prot_length = item
            accession = ida[0].split(':')[0]
            domains_to_check[accession] = accessions[accession]
            domains_to_check[accession]["domains"] = ida[0]
            domains_to_check[accession]["count_ida"] = count_ida
            domains_to_check[accession]["protein_count"] = count_total
            domains_to_check[accession]["representative"] = f"{protein}:{dom_length}/{prot_length}"
    
    # print(domains_to_check)
    print("pfam_acc,pfam_type,ipr_acc,ipr_type,most_common_ida,protein_with_ida,total_protein_count,representative_protein:dom_length/prot_length")
    for accession, info in domains_to_check.items():
        print(f"{accession},{info['m_type']},{info['e_acc']},{info['e_type']},{info['domains']},{info['count_ida']},{info['protein_count']},{info['representative']}")