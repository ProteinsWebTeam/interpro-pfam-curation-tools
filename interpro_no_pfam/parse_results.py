from search_entry_no_pfam import getOracleConnection
import argparse
import os, re
from configparser import ConfigParser
import cx_Oracle
import traceback
import subprocess
import shutil

def get_pmid(pub_id, cursor):
    request = f"""SELECT pubmed_id FROM CITATION WHERE pub_id=:1"""
    cursor.execute(request, (pub_id,))
    pubmed_id = cursor.fetchone()
    if pubmed_id is None:
        return None
    else:
        return pubmed_id[0]
    
def get_uniprot_ids(entry, cursor, outputdir):

    output_file = os.path.join(outputdir, f"{entry}.tsv")
    if os.path.exists(output_file):
        return
    
    request = f"""SELECT m.PROTEIN_AC, m.pos_from, m.pos_to
                    FROM MATCH m
                    JOIN ENTRY2METHOD e ON m.METHOD_AC = e.METHOD_AC
                    WHERE e.ENTRY_AC=:1
                    order by m.protein_ac
                """
    cursor.execute(request, (entry,))
    matches = cursor.fetchall()

    save_protein = None
    save_start = None
    save_end = None

    with open(output_file, 'w') as outfile:
        for match in matches:
            print(match)
            if save_protein is None:
                save_protein = match[0]
                save_start = match[1]
                save_end = match[2]
            else:
                if save_protein == match[0]:
                    if match[1] < save_start:
                        save_start = match[1]
                    elif match[2] > save_end:
                        save_end = match[2]
                else:
                    outfile.write(f"{save_protein}\t{save_start}\t{save_end}\n")
                    save_protein = match[0]
                    save_start = match[1]
                    save_end = match[2]

        outfile.write(f"{save_protein}\t{save_start}\t{save_end}\n")


def overlap(start1, end1, start2, end2):
    """Check if two ranges overlap.""" 
    return (start1 <= start2 and end1 <= end2 and start2 <= end1) or (start2 <= start1 and end2 <= end1 and start1 <= end2) or (start2 <= start1 and end2 >= end1) or (start2 >= start1 and end2 <= end1)

import subprocess

def replace_greek_letters(text):
    # Mapping of Greek letters to words
    greek_map = {
        'α': 'alpha',
        'β': 'beta',
        'γ': 'gamma',
        'δ': 'delta',
        'ε': 'epsilon',
    }
    
    greek_regex = re.compile('|'.join(re.escape(key) for key in greek_map.keys()))
    
    # Replace each Greek letter or character with its corresponding word
    return greek_regex.sub(lambda match: greek_map[match.group(0)], text)


def clean_html(raw_html, cursor, citation_counter, pmid_list):
    # Remove <p> tags and replace with newline
    raw_html = re.sub(r'<p>\s*', '', raw_html)
    raw_html = re.sub(r'</p>\s*', '\n', raw_html)

    # Remove <ul> tags
    raw_html = re.sub(r'<ul>\s*', '', raw_html)
    raw_html = re.sub(r'</ul>\s*', '', raw_html)

    # Replace <li> tags and handle commas
    def replace_li(match):
        content = match.group(1).strip()
        if content.endswith('.'):
            return content
        return f"{content}, "

    raw_html = re.sub(r'<li>(.*?)</li>', replace_li, raw_html, flags=re.DOTALL)

    # Remove other HTML tags
    cleaned_text = re.sub(r'<.*?>', '', raw_html)
    cleaned_text = replace_greek_letters(cleaned_text)

    # Remove blank lines
    cleaned_text = "\n".join(line for line in cleaned_text.splitlines() if line.strip())

    # Function to replace each [[cite:PUB]] with a reference number
    def replace_citations(match):
        nonlocal citation_counter
        nonlocal pmid_list
        citations = match.group(0)

        def replace_single_citation(pub_match):
            nonlocal citation_counter
            nonlocal pmid_list
            pub_id = pub_match.group(1)
            pmid = get_pmid(pub_id, cursor)

            # Replace citation with number
            if pmid in pmid_list:
                ref_number = pmid_list[pmid]
            else:
                ref_number = citation_counter
                citation_counter += 1
                pmid_list[pmid] = ref_number
            return str(ref_number)
        
        # Replace each PUB within this group
        return re.sub(r'\[cite:(.*?)\]', replace_single_citation, citations)

    # Replace consecutive citations
    cleaned_text = re.sub(r'\[\[cite:.*?\](?:,\s*\[\[cite:.*?\])*\]', replace_citations, cleaned_text)

    return cleaned_text, citation_counter, pmid_list


def format_annotation(annotation):
    lines = []
    current_line = "CC   "
    
    for word in annotation.split():
        if len(current_line) + len(word) + 1 > 80:  # +1 for space
            lines.append(current_line)
            current_line = "CC   " + word
        else:
            if current_line != "CC   ":
                current_line += " "
            current_line += word

    lines.append(current_line)  # Add last line
    return "\n".join(lines)


def get_interpro_data(interpro_ac, ora_cursor, outputdir):
    output_file = os.path.join(outputdir, f"{interpro_ac}.desc")
    if os.path.exists(output_file):
        return

    request = f"""SELECT e.name, e.short_name, e2c.order_in, c.text
                    FROM ENTRY e
                    JOIN ENTRY2COMMON e2c ON e.ENTRY_AC = e2c.ENTRY_AC
                    JOIN COMMON_ANNOTATION c ON e2c.ANN_ID = c.ANN_ID
                    WHERE e.ENTRY_AC=:1
                    ORDER BY e2c.order_in
                """
    ora_cursor.execute(request, (interpro_ac,))
    name = ""
    short_name = ""
    annotation = ""
    ref_counter = 1
    pmid_list = {}
    
    for row in ora_cursor.fetchall():
        if name == "":
            name = f"DE   {row[0]}"
            short_name = f"ID   {row[1]}"
        
        annot, ref_counter, pmid_list = clean_html(row[3], ora_cursor, ref_counter, pmid_list)
        annotation += annot + " "
    # print(annotation)
    formatted_annotation = format_annotation(annotation.strip())

    # print(formatted_annotation)
    desc_file = os.path.join(outputdir, "DESC")
    with open (desc_file, 'w') as outfile:        
        outfile.write(f"{short_name}\n{name}\n")

    # print(pmid_list)

    os.chdir(outputdir)
    return_code = None
    with open(os.path.join(outputdir,f"add_ref.log"),"wb") as f:
        for pmid in pmid_list.keys():
            cmd = ["perl", "/hps/software/users/agb/pfam/software/bin/add_ref", str(pmid)]
            popen = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, shell=False)
            return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    if os.path.exists(desc_file):
        with open(desc_file, "a") as f:
            f.write(formatted_annotation)

        shutil.move(desc_file, output_file)
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    oracle_connectionString = config["database"]["ippro-interpro"]
    connection, cursor = getOracleConnection(oracle_connectionString)

    results = {"count60+":0, "count50-59":0, "count40-49":0, "count30-39":0, "count20-29":0, "count10-19":0, "count1-9":0, "count0":0}

    accessions0dir = config["files"]["output_uniprot_0overlap"]
    accessions10dir = config["files"]["output_uniprot_10overlap"]
    output_allfile = config["files"]["outputfile"]

    os.makedirs(accessions0dir, exist_ok=True)
    os.makedirs(accessions10dir, exist_ok=True)

    # get_interpro_data("IPR005746", cursor, accessions0dir)

    with open(output_allfile, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            entry, nb_proteins, nb_proteins_in_pfam, nb_pfams = line.split('\t')

            percentage = int(nb_proteins_in_pfam)*100/int(nb_proteins)
            print(entry, percentage)
            if percentage == 0:
                results['count0'] += 1
                get_uniprot_ids(entry, cursor, accessions0dir)
                get_interpro_data(entry, cursor, accessions0dir)
            elif percentage < 10:
                results['count1-9'] += 1
                get_uniprot_ids(entry, cursor, accessions10dir)
            elif percentage < 20:
                results['count10-19'] += 1
            elif percentage < 30:
                results['count20-29'] += 1
            elif percentage < 40:
                results['count30-39'] += 1
            elif percentage < 50:
                results['count40-49'] += 1
            elif percentage < 60:
                results['count50-59'] += 1
            else:
                results['count60+'] += 1
            
    cursor.close()
    connection.close()

    print(results)