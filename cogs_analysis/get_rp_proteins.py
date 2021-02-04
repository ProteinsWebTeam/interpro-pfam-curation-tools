import redis


def redis_connection():
    # Connection to REDIS
    REDIS_HOST = "ebi-pdbe.ebi.ac.uk"
    REDIS_PORT = "54321"

    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT)


server = redis_connection()
redis_rp = 'jaccard_cog_rp'

server.delete(redis_rp)

# load representative proteins into redis queue
print("Loading representative proteins into memory")
with open("reference_proteome_accessions_2.txt", 'r') as RP:
    for line in RP:
        protein = line.strip('\n')
        server.hset(redis_rp, protein, 1)
