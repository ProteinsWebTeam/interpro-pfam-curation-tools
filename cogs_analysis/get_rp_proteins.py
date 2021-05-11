def getProteinAccessions(server, ref_prot_file, redis_rp):
    # load representative proteins into redis queue
    print("Loading representative proteins into memory")
    with open(ref_prot_file, "r") as RP:
        for line in RP:
            protein = line.strip("\n")
            server.hset(redis_rp, protein, 1)

