import cx_Oracle
import redis
import os
import subprocess
import json
from shutil import rmtree
from itertools import groupby


def redis_connection():
    # Connection to REDIS
    REDIS_HOST = "ebi-pdbe.ebi.ac.uk"
    REDIS_PORT = "54321"

    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT)


def getMidpoints(fragmentlist):
    fragments = fragmentlist.split(',')
    midpoints = []

    for frag in fragments:
        start, end = frag.split('-')
        midpoints.append((int(end)-int(start))/2+float(start))

    return midpoints


def main():
    # redis connection
    server = redis_connection()
    redis_protein = 'jaccard_cogs_prot'
    redis_rp = 'jaccard_cog_rp'
    redis_ipr_count = 'jaccard_cog_ipr'

    server.delete(redis_protein)
    server.delete(redis_ipr_count)

    # iprel connection
    dsn = cx_Oracle.makedsn("***REMOVED***", "***REMOVED***", 'iprel')
    connection = cx_Oracle.connect(
        user='interpro_webserver', password='iwebserv', dsn=dsn)
    dbcursor = connection.cursor()

    # get ipr and representative proteins
    print("Calculating InterPro integrated entries midpoints")

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
                    ipr_count = int(server.hget(
                        redis_ipr_count, ipr).decode('utf-8'))
                    server.hset(redis_ipr_count, ipr, ipr_count+1)
                else:
                    server.hset(redis_ipr_count, ipr, 1)

                # add ipr midpoint to protein dict
                for midpoint in midpoints:
                    protein_list.setdefault(ipr, []).append(midpoint)

            server.hset(redis_protein, group[0], json.dumps(protein_list))

    print("Calculating InterPro integrated entries midpoints completed")


main()
