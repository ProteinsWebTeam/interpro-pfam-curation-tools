import os
from shutil import rmtree

import cx_Oracle
import redis


def cleanDir(outputdir):
    # clean result directories
    try:
        rmtree(outputdir)
    except:
        pass

    os.makedirs(outputdir, exist_ok=True)


def redisConnection(config):
    # Connection to REDIS
    REDIS_HOST = config["redis"]["redis_host"]
    REDIS_PORT = config["redis"]["redis_port"]

    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT)


def dbConnection(config, schema):

    host = config[schema]["host"]
    port = config[schema]["port"]
    user = config[schema]["user"]
    password = config[schema]["password"]
    print(user, password, host)
    dsn = cx_Oracle.makedsn(host, port, schema)

    connection = cx_Oracle.connect(user=user, password=password, dsn=dsn)
    dbcursor = connection.cursor()

    return connection, dbcursor


def closeConnection(connection, cursor):
    cursor.close()
    connection.close()
