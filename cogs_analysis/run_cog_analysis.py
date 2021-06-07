import argparse
import os, sys
from configparser import ConfigParser

from utils import dbConnection, redisConnection, closeConnection, cleanDir
from get_cog_list import getCOGlist, getOverlap
from get_rp_proteins import getProteinAccessions
from IPR_RP_AC import getRepresentativeInfo

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    parser.add_argument(
        "-1",
        "--step1",
        action="store_true",
        help="load list of representative proteins into memory",
    )
    parser.add_argument(
        "-2", "--step2", action="store_true", help="process list of representative proteins"
    )

    parser.add_argument("-3", "--step3", action="store_true", help="load list of COGs into memory")
    parser.add_argument(
        "-4",
        "--step4",
        action="store_true",
        help="search intersect between COGs and representative proteins",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}'': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    pfile = config["files"]["pfile"]
    if not os.path.isfile(pfile):
        parser.error(f"cannot open {pfile}: " f"no such file or directory")

    server = redisConnection(config)
    redis_rp = config["redis"]["redis_rp"]
    redis_protein = config["redis"]["redis_protein"]
    redis_ipr_count = config["redis"]["redis_ipr_count"]

    if args.step1:
        print("Starting step 1")
        print(f"Deleting previous data from redis queue {redis_rp}")
        server.delete(redis_rp)

        getProteinAccessions(server, pfile, redis_rp)

    if args.step2:
        if server.exists(redis_rp):
            print("Starting step 2")

            # iprel connection
            ##### This need to be replaced as IPREL schema and the table used don't exist anymore ######
            connection, dbcursor = dbConnection(config, "iprel")

            getRepresentativeInfo(config, server, dbcursor)
            closeConnection(connection, dbcursor)
        else:
            print(f"Error redis queue '{redis_rp}' empty, restart from Step 1")
            sys.exit(1)

    if args.step3:
        print("Starting step 3")
        # iptst connection
        connection, dbcursor = dbConnection(config, "iptst")

        redis_cog = config["redis"]["jaccard_cogs"]
        print(f"Deleting previous data from redis queue {redis_cog}")
        server.delete(redis_cog)

        outputdir_log = config["files"]["outputdir_log"]
        outputdir_files = config["files"]["outputdir_files"]
        cleanDir(outputdir_log)
        cleanDir(outputdir_files)

        getCOGlist(config, server, dbcursor)
        closeConnection(connection, dbcursor)
    if args.step4:
        if server.exists(redis_protein):
            print("Starting step 4")
            # getOverlap(outputdir_log, outputdir_files)
        else:
            if server.exists(redis_rp):
                print(f"Error redis queue '{redis_protein}' empty, restart from Step 2")
            else:
                print(f"Error redis queue '{redis_protein}' empty, restart from Step 1")
            sys.exit(1)
