import cx_Oracle
import redis
import argparse
import json


def redis_connection():
    # Connection to REDIS
    REDIS_HOST = "ebi-pdbe.ebi.ac.uk"
    REDIS_PORT = "54321"

    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT)


def isIntersect(midpoint, start, end):
    if midpoint >= start and midpoint <= end:
        return 1
    else:
        return 0


def main():
    # redis connection
    server = redis_connection()

    # iptst connection
    dsn = cx_Oracle.makedsn("ora-dlvm5-014.ebi.ac.uk", "1521", 'iptst')
    connection_cog = cx_Oracle.connect(
        user='OPS$typhaine', password='typhaine55', dsn=dsn)
    # connection_cog = cx_Oracle.connect('ops$hychang/changeme@iptst')
    dbcursor_cog = connection_cog.cursor()

    cog = 1

    while cog:
        # get IPR/signature from the REDIS_QUEUE
        cog = server.lpop(ARGS.redis_cog)

        if not cog:
            break

        cog = str(cog.decode('utf8'))

        # get common protein accessions and intersect between COGs and IPR
        print("Searching for intersecting proteins betwwen COGS and InterPro")

        request_cog = "SELECT PROTEIN_AC,METHOD_AC,POS_FROM,POS_TO FROM INTER_WORK.MATCH_COG_20171211_ALL where method_ac=:1"
        dbcursor_cog.execute(request_cog, (cog,))
        member1list = dbcursor_cog.fetchall()

        intersect_dict = {}
        cogs_count = {}

        for row in member1list:
            protein, cog, start, end = row
            if server.hexists(ARGS.rp_proteins, protein):
                # initialize protein count for cogs
                if cog in cogs_count:
                    cogs_count[cog].add(protein)
                else:
                    cogs_count[cog] = set([protein])

                # if the protein in found with at least one IPR
                if server.hexists(ARGS.protein_list, protein):
                    list_ipr = json.loads(server.hget(
                        ARGS.protein_list, protein).decode('utf8'))

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
        with open(ARGS.output, 'a') as file_to_print:

            for key in sorted(intersect_dict.keys()):
                intersect = len(intersect_dict[key])
                (cog, ipr) = key.split('_')

                set1 = len(cogs_count[cog])
                set2 = int(server.hget(
                    ARGS.ipr_count, ipr).decode('utf-8'))

                union = set1+set2-intersect

                ji = round(intersect/union, 2)
                jc1 = round(intersect/float(set1), 2)
                jc2 = round(intersect/float(set2), 2)

                file_to_print.write("{},{},{},{},{},{},{},{},{}\n".format(
                    cog, ipr, ji, jc1, jc2, intersect, union, set1, set2))

    print("Calculation completed")


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--redis_cog", help='name of the redis queue for cogs accessions', default='jaccard_cogs')
    parser.add_argument(
        "--rp_proteins", help='name of the redis queue for representative proteins', default='jaccard_cog_rp')
    parser.add_argument(
        "--protein_list", help='name of the redis queue for ipr and midpoints per protein', default='jaccard_cogs_prot')
    parser.add_argument(
        "--ipr_count", help='name of the redis queue for total number of proteins per ipr', default='jaccard_cog_ipr')
    parser.add_argument(
        "output", help='name of the output file for the current job')

    ARGS = parser.parse_args()

    main()
