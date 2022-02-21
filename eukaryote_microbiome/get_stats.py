import re, os
import argparse
from multiprocessing import Pool
from functools import partial
import itertools


def process(cluster_dir, cluster):
    countswiss = 0
    count_total = 0
    content = ""
    rep = ""
    counting = ""

    for item in cluster:
        if re.match(r"^>", item):
            if re.match(r">sp", item):  # seq from SwissProt
                countswiss += 1
                acc = item.split("|")[1]
                item = item.replace("sp|", "")
            else:  # seq from TrEMBL
                acc = item.split("|")[1]
                item = item.replace("tr|", "")

            # total number of seq
            count_total += 1

            # update rep to new representative accession
            if count_total == 1:
                rep = acc
            item = item.replace("|", " ")
        content += item

    # for clusters with at least 2 sequences
    if count_total > 1:
        # generate % SwissProt sequences found in cluster
        percentswiss = round(countswiss * 100 / count_total, 2)
        subdir = os.path.join(cluster_dir, rep[0:3])
        os.makedirs(subdir, exist_ok=True)

        # save clusters in different files
        with open(os.path.join(subdir, f"{rep}.fa"), "w") as clusterf:
            clusterf.write(content)

        counting = f"{rep}\t{count_total}\t{percentswiss}\n"
        return counting


def cluster_separator(line):
    if re.match(r"^>[A-Z0-9]+$", line):
        return line


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputfile", help="file countaining cluster_sequences (fasta)", required=True
    )

    args = parser.parse_args()

    cluster_dir = os.path.join(os.path.dirname(args.inputfile), "clusters")

    os.makedirs(cluster_dir, exist_ok=True)

    print("Getting clusters' statistics")

    num_chunks = 5
    results = []
    counter = 0

    if os.path.isfile(args.inputfile):
        outputfile = os.path.join(cluster_dir, f"{args.inputfile}_percent_euk")
        with open(args.inputfile, "r") as inputf, open(outputfile, "w") as output:
            chunks = itertools.groupby(inputf, cluster_separator)
            try:
                while True:
                    groups = []
                    for key, chunk in itertools.islice(chunks, num_chunks):
                        if not key:
                            groups.append(list(chunk))
                    if groups:
                        with Pool(5) as pool:
                            for cluster in pool.map(partial(process, cluster_dir), groups):
                                if cluster != None:
                                    counter += 1
                                    output.write(cluster)
                                # end search if 10K clusters found
                                if counter >= 10000:
                                    raise ValueError("Stop!!!")
                    else:
                        break
            except ValueError:
                print(f"Clustering check done, {counter} clusters saved")

