import json
from itertools import groupby


def getMidpoints(fragmentlist):
    fragments = fragmentlist.split(",")
    midpoints = []

    for frag in fragments:
        start, end = frag.split("-")
        midpoints.append((int(end) - int(start)) / 2 + float(start))

    return midpoints


def getRepresentativeInfo(config, server, dbcursor):

    redis_protein = config["redis"]["redis_protein"]
    redis_rp = config["redis"]["redis_rp"]
    redis_ipr_count = config["redis"]["redis_ipr_count"]

    server.delete(redis_protein)
    server.delete(redis_ipr_count)

    # get ipr and representative proteins
    print("Calculating InterPro integrated entries midpoints")

    ##### This need to be replaced as IPREL schema and this table don't exist anymore ######
    #### idea: use InterPro API: https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/ #####
    request = "SELECT s.ENTRY_AC,s.PROTEIN_AC,s.FRAGMENTS \
        FROM INTERPRO.UPI_SUPERMATCH_STG_NO_UPI_NEW s \
        JOIN INTERPRO.ENTRY e on e.ENTRY_AC=s.ENTRY_AC \
        WHERE (e.ENTRY_TYPE='F' or e.ENTRY_TYPE='D') and e.CHECKED='Y' \
        ORDER BY s.protein_ac,s.entry_ac"

    dbcursor.execute(request)

    ipr_count = {}

    for group in groupby(dbcursor, lambda x: x[1]):
        # verify if protein in representative list
        if server.hexists(redis_rp, group[0]):
            protein_list = {}
            for row in group[1]:
                ipr = row[0]
                midpoints = getMidpoints(row[2])

                # add protein to ipr_count list
                if server.hexists(redis_ipr_count, ipr):
                    ipr_count = int(server.hget(redis_ipr_count, ipr).decode("utf-8"))
                    server.hset(redis_ipr_count, ipr, ipr_count + 1)
                else:
                    server.hset(redis_ipr_count, ipr, 1)

                # add ipr midpoint to protein dict
                for midpoint in midpoints:
                    protein_list.setdefault(ipr, []).append(midpoint)

            server.hset(redis_protein, group[0], json.dumps(protein_list))

    print("Calculating InterPro integrated entries midpoints completed")

