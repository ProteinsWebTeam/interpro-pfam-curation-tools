import cx_Oracle
import redis
import os
import subprocess
from shutil import rmtree


def redis_connection():
    # Connection to REDIS
    REDIS_HOST = "ebi-pdbe.ebi.ac.uk"
    REDIS_PORT = "54321"

    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT)


def main():
    # clean result directories
    outputdir_log = "/hps/nobackup/production/interpro/jaccard_cogs/log"
    outputdir_files = "/hps/nobackup/production/interpro/jaccard_cogs/files"

    try:
        rmtree(outputdir_log)
        rmtree(outputdir_files)
    except:
        pass

    os.makedirs(outputdir_log, exist_ok=True)
    os.makedirs(outputdir_files, exist_ok=True)

    server = redis_connection()
    redis_cog = 'jaccard_cogs'
    server.delete(redis_cog)

    dsn = cx_Oracle.makedsn("ora-dlvm5-014.ebi.ac.uk", "1521", 'iptst')
    connection_cog = cx_Oracle.connect(
        user='OPS$typhaine', password='typhaine55', dsn=dsn)
    # connection_cog = cx_Oracle.connect('ops$hychang/changeme@iptst')
    dbcursor_cog = connection_cog.cursor()

    # get cog list from db
    print("Searching for COGS accessions")

    request_cog = "SELECT distinct METHOD_AC FROM INTER_WORK.MATCH_COG_20171211_ALL"
    # WHERE METHOD_AC='COG0409'
    dbcursor_cog.execute(request_cog)
    member1list = dbcursor_cog.fetchall()

    for row in member1list:
        server.lpush(redis_cog, row[0])

    dbcursor_cog.close()
    connection_cog.close()

    # overlap
    script_path = os.path.dirname(os.path.realpath(__file__))
    overlapscript = os.path.join(script_path, 'search_intersect.py')

    print("Search overlapping entries")

    log_file = os.path.join(outputdir_log, "overlap_%I.log")

    subprocess.run(['bsub', '-o', log_file, '-M', '8192', '-R', '"rusage[mem=8192]"', '-J',
                    'overlap_cogs[1-8]', '/nfs/msd/sw/RedHatEnterpriseServer-7.3-x86_64/bin/python3', overlapscript, os.path.join(outputdir_files, 'overlap$LSB_JOBINDEX')])


main()
