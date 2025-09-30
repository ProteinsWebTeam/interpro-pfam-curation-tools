import mysql.connector
import csv
import os
import re

def sql_connection(dbinfo):
    connection = mysql.connector.connect(
        host=dbinfo[0],
        user=dbinfo[1],
        password=dbinfo[2],
        database=dbinfo[3],
        port=dbinfo[4]
    )

    # Check if the connection is successful
    if connection.is_connected():
        pass
    else:
        print("Failed to connect to the database")
        exit()
    
    return connection


def write_results_to_tsv(file_name, pfam_data, results):
    """Write the results to a TSV file with the specified columns."""
    with open(file_name, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')

        # Write the header
        # writer.writerow(['initial_pfam_acc', 'initial_prot','seq_version', 'start', 'end', 'representative', 'new_pfam_acc', 'new_clan', 'new_prot', 'start', 'end'])

        try:
            pfam_acc, protein_acc, seq_version, location = pfam_data.split("_")
        except ValueError:
            pfam_acc, location = pfam_data.split("_")
            protein_acc = None
            seq_version = None
        start, end = location.split("-")
        start = start.replace("res", "")

        for row in results:
            to_write = [pfam_acc, protein_acc, seq_version, start, end]

            for value in row.values():
                to_write.append(value)

            writer.writerow(to_write)


def extract_pfam_representatives_from_file(chopped_struct_dir, pfam_id=None):
    """Extract PF accession numbers and protein identifiers from files in a given directory as a list of dictionaries."""
    results = {}
    pattern = r'^(PF\d{5})_(\w+)_(\d{1})_res(\d+)-(\d+)\.cif$'  # Updated regex pattern to also capture seq_start and seq_end

    outfile = None
    if os.path.isdir(chopped_struct_dir):
        for filename in os.listdir(chopped_struct_dir):
            if os.path.isfile(os.path.join(chopped_struct_dir, filename)):
                match = re.match(pattern, filename)
                if match:
                    pfam_acc = match.group(1)
                    pfamseq_acc = match.group(2)
                    seq_version = int(match.group(3))
                    seq_start = int(match.group(3))
                    seq_end = int(match.group(4))

                    results[pfam_acc] = {
                        'pfamA_acc': pfam_acc,
                        'pfamseq_acc': pfamseq_acc,
                        'seq_version': seq_version,
                        'seq_start': seq_start,
                        'seq_end': seq_end
                    }
                    if pfam_id == pfam_acc:
                        outfile = os.path.join(chopped_struct_dir, filename)

    return results, outfile