import os
import subprocess


def getCOGlist(redis_cog, server, dbcursor):
    # get cog list from db
    print("Searching for COGS accessions")

    request_cog = "SELECT distinct METHOD_AC FROM INTER_WORK.MATCH_COG_20171211_ALL"
    # WHERE METHOD_AC='COG0409'
    dbcursor.execute(request_cog)
    member1list = dbcursor.fetchall()

    for row in member1list:
        server.lpush(redis_cog, row[0])


def getOverlap(outputdir_log, outputdir_files):
    # overlap
    script_path = os.path.dirname(os.path.realpath(__file__))
    overlapscript = os.path.join(script_path, "search_intersect.py")

    print("Search overlapping entries")

    log_file = os.path.join(outputdir_log, "overlap_%I.log")

    subprocess.run(
        [
            "bsub",
            "-o",
            log_file,
            "-M",
            "8192",
            "-R",
            '"rusage[mem=8192]"',
            "-J",
            "overlap_cogs[1-8]",
            "/nfs/msd/sw/RedHatEnterpriseServer-7.3-x86_64/bin/python3",
            overlapscript,
            os.path.join(outputdir_files, "overlap$LSB_JOBINDEX"),
        ]
    )

