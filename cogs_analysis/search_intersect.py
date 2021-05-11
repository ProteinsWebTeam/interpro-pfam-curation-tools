import argparse
import json
import os
from configparser import ConfigParser

from cogs_analysis.utils import closeConnection, dbConnection, redisConnection


def isIntersect(midpoint, start, end):
    if midpoint >= start and midpoint <= end:
        return 1
    else:
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file", required=True)
    parser.add_argument("output", help="name of the output file for the current job")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    args = parser.parse_args()

    # redis connection
    server = redisConnection(config)
    redis_cog = config["redis"]["redis_cog"]
    rp_proteins = config["redis"]["redis_rp"]
    protein_list = config["redis"]["redis_protein"]
    ipr_count = config["redis"]["redis_ipr_count"]

    # iptst connection
    connection_cog, dbcursor_cog = dbConnection(config, "iptst")

    cog = 1

    while cog:
        # get IPR/signature from the REDIS_QUEUE
        cog = server.lpop(redis_cog)

        if not cog:
            break

        cog = str(cog.decode("utf8"))

        # get common protein accessions and intersect between COGs and IPR
        print("Searching for intersecting proteins betwwen COGS and InterPro")

        request_cog = "SELECT PROTEIN_AC,METHOD_AC,POS_FROM,POS_TO FROM INTER_WORK.MATCH_COG_20171211_ALL where method_ac=:1"
        dbcursor_cog.execute(request_cog, (cog,))
        member1list = dbcursor_cog.fetchall()

        intersect_dict = {}
        cogs_count = {}

        for row in member1list:
            protein, cog, start, end = row
            if server.hexists(rp_proteins, protein):
                # initialize protein count for cogs
                if cog in cogs_count:
                    cogs_count[cog].add(protein)
                else:
                    cogs_count[cog] = set([protein])

                # if the protein in found with at least one IPR
                if server.hexists(protein_list, protein):
                    list_ipr = json.loads(server.hget(protein_list, protein).decode("utf8"))

                    # for each IPR matching the protein, searching if they intersect with the current COG
                    for ipr in list_ipr:
                        for midpoint in list_ipr[ipr]:
                            if isIntersect(float(midpoint), start, end):
                                key = "{}_{}".format(cog, ipr)

                                if key in intersect_dict:
                                    intersect_dict[key].add(protein)
                                else:
                                    intersect_dict[key] = set([protein])

                                break

        print("Saving results in file")
        with open(args.output, "a") as file_to_print:

            for key in sorted(intersect_dict.keys()):
                intersect = len(intersect_dict[key])
                (cog, ipr) = key.split("_")

                set1 = len(cogs_count[cog])
                set2 = int(server.hget(ipr_count, ipr).decode("utf-8"))

                union = set1 + set2 - intersect

                ji = round(intersect / union, 2)
                jc1 = round(intersect / float(set1), 2)
                jc2 = round(intersect / float(set2), 2)

                file_to_print.write(
                    "{},{},{},{},{},{},{},{},{}\n".format(
                        cog, ipr, ji, jc1, jc2, intersect, union, set1, set2
                    )
                )

    print("Calculation completed")
    closeConnection(connection_cog, dbcursor_cog)

