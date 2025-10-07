from configparser import ConfigParser
import argparse
import os
import pandas as pd
import oracledb
from datetime import datetime

# function to get the information from the InterPro database

def get_updated_DUFs (ip_user, ip_pw, ip_db, ip_pf, ip_lf):
    ip_query = f"""
    select n.ip_acc, o.old_short_name, o.old_name, n.new_short_name, n.new_name, n.last_edit, c.pubmed_id as pmid, p.protein_ac as SwissProt
    from 
    (
        select *
        from
            (
            select ea.entry_ac, max(ea.short_name) as old_short_name, max(ea.name) as old_name, max(ea.timestamp) as previous_edit
            from interpro.entry_audit ea
            where ea.checked='Y'
            and cast(ea.timestamp as date) < '{ip_pf}'
            group by ea.entry_ac
            ) eam
        where (eam.old_short_name LIKE '%DUF%'
        or eam.old_name LIKE '%UPF%'
        or eam.old_name LIKE '%UCP%')
    ) o
    right join
        (
        select e.entry_ac as IP_acc, e.name as new_name, e.short_name as new_short_name, cast(e.timestamp as date) as last_edit
        from interpro.entry e
        where e.checked='Y'
        and cast(e.timestamp as date) >= '{ip_pf}'
        ) n on o.entry_ac=n.ip_acc
            join interpro.entry2pub e2p on e2p.entry_ac=n.ip_acc
    join INTERPRO.citation c on c.pub_id=e2p.pub_id
    join INTERPRO.entry2method e2m on e2m.entry_ac=n.ip_acc
    join INTERPRO.match m on m.method_ac=e2m.method_ac
    join INTERPRO.protein p on p.protein_ac=m.protein_ac
    where o.entry_ac is not null
    and cast(o.previous_edit as date) < '{ip_lf}'
    and n.new_short_name not LIKE '%DUF%'
    and p.dbcode='S'

    """

    try:
        ip_connection = oracledb.connect(
            user=ip_user,
            password=ip_pw,
            dsn=ip_db)
            
        try:
            with open('updated_DUFs.log', 'a') as f:
                f.write(f"Successfully connected to the InterPro Database. \n \n")
        except:
            pass

        ip_df = pd.read_sql_query(ip_query, ip_connection)

        grouping_cols = ['IP_ACC', 'OLD_SHORT_NAME','OLD_NAME','NEW_SHORT_NAME','NEW_NAME','LAST_EDIT']
        agg_cols = ['PMID', 'SWISSPROT']

        ip_df = ip_df.groupby(grouping_cols)[agg_cols].agg(lambda x: x.dropna().unique().tolist()).reset_index()

        ip_connection.close()

    except:
        print("Couldn't connect to the InterPro database.")
        ip_df = pd.DataFrame()
        try:
            with open('updated_DUFs.log', 'a') as f:
                f.write(f"Couldn't connect to the InterPro database. \n \n")
        except:
            pass


    try:
        with open('updated_DUFs.log', 'a') as f:
            f.write(f"ipdf:\n {ip_df.head(3)} \n \n")
    except:
        pass

    return ip_df


###
# Beatriz Lazaro Pinto (June 2025)
#
# This script gets a csv file with a table of old InterPro DUFs that have been renamed from the last InterPro release 
# with the PMIDs and the Swissprot proteins currently linked to the entry 
# 
# The date of the previous and last Pronto freeze should be specified in the config_last_updated_InterPro_DUFs.ini file
# together with the parameters of the database
#
# The config_last_updated_InterPro_DUFs.ini file must be in the folder where this script is run
###

if __name__ == "__main__":

    pd.set_option('display.max_columns', None)

    # create log file
    with open('updated_DUFs.log', 'w') as f:
        start_time = datetime.now()
        f.write(f'Starting time: {str(start_time)} \n \n')

    # get information from config file
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f'Cannot open {args.config}: no such file or directory')

    config = ConfigParser()
    config.read(args.config)

    # get ip info and create df of the DUFs in IP with a Pfam contributing signature

    ip_user = str(config['interpro']['user'])
    ip_pw = str(config['interpro']['password'])
    ip_db = str(config['interpro']['db'])

    ip_pf = str(config['release']['date_previous_freeze'])
    ip_lf = str(config['release']['date_last_freeze'])
    
    updates = get_updated_DUFs(ip_user, ip_pw, ip_db, ip_pf, ip_lf)

    updates.to_csv('DUF_updates.csv')
    print("File succesfully created.")

    with open('updated_DUFs.log', 'a') as f:
        f.write(f"File succesfully created. \n \n")
        # record end time and total time
        final_time = datetime.now() - start_time
        f.write(f'Total time employed: {final_time} minutes.')