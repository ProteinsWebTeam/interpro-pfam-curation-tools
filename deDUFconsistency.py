from configparser import ConfigParser
import argparse
import os
import pandas as pd
import oracledb
import mysql.connector

pd.set_option('display.max_columns', None)

# function to get the list of DUFs in pfam (from df of all Pfam families)
def get_pf_duf_df(pf_df):
    pf_duf_df = pf_df[pf_df['pfamA_id'].str.contains('DUF')]
    return pf_duf_df

# function to get a list of Pfam accession that are not in a list of DUFs in InterPro
def check_mismatching(ip_df, pf_duf_df):
    list_mismatches = []

    list_pfam_duf_ac = pf_duf_df.index.values.tolist()
    print(len(list_pfam_duf_ac), list_pfam_duf_ac[:5])
    list_ip_duf_ac = ip_df['METHOD_AC'].tolist()
    print(len(list_ip_duf_ac), list_ip_duf_ac[:5])

    list_mismatches = list(set(list_ip_duf_ac).symmetric_difference(set(list_pfam_duf_ac)))
    print('list1: ',len(list_mismatches), list_mismatches[:5])
    
    list_mismatches = [i for i in list_mismatches if i[2:] <= '23022']
    print('list2: ',len(list_mismatches), list_mismatches[:5])

    return list_mismatches

# create dataframe from the list of mismatches with all the information
# that can help quickly see if further curation is worth it
def create_df(list_mismatches, ip_df, pf_df):

    # create dataframe with InterPro DUFs that are not DUFs in Pfam
    df_ip_dufs = pd.DataFrame(columns=('IP_ac', 'IP_short_name', 'IP_name', 'Pfam_ac', 'Pfam_ID', \
                    'Pfam_description', 'Pfam_previous_ID'))
    df_ip_dufs.Pfam_ac = list_mismatches

    list_pf_dufs_with_diff_ip_name =[]

    for i, row in df_ip_dufs.iterrows():
        # exclude cases where the entry is a ~DUF in Pfam but not in IP (so it's not in ip_df)
        try:
            n_row_ip = int(list(ip_df.index[ip_df['METHOD_AC'] == str(list_mismatches[i])])[0])
            df_ip_dufs['IP_ac'][i] = ip_df['ENTRY_AC'][n_row_ip]
            df_ip_dufs['IP_short_name'][i] = ip_df['SHORT_NAME'][n_row_ip]
            df_ip_dufs['IP_name'][i] = ip_df['IP_NAME'][n_row_ip]
            df_ip_dufs['Pfam_ID'][i] = ip_df['PF_ID'][n_row_ip]
            df_ip_dufs['Pfam_description'][i] = ip_df['DESCRIPTION'][n_row_ip]
            df_ip_dufs['Pfam_previous_ID'][i] = pf_df['previous_id'][list_mismatches[i]]
        except:
            list_pf_dufs_with_diff_ip_name.append(list_mismatches[i])

    df_ip_dufs = df_ip_dufs.dropna(axis=0,thresh=2)

    # create dataframe with Pfam DUFs that are not DUFs in InterPro

    df_pf_dufs = pd.DataFrame()

    for i in list_pf_dufs_with_diff_ip_name[:5]:
        ip_query_list = f"""
        select e.entry_ac, e2m.method_ac, e.name as IP_NAME, e.short_name, em.description, em.name as PF_ID
        from interpro.entry e
        join interpro.entry2method e2m on e.entry_ac=e2m.entry_ac
        join interpro.method em on e2m.method_ac=em.method_ac
        where e.checked='Y'
        and e2m.method_ac = '{i}'
        """
        list_pf_could_not_get_data = []
        try:
            ip_new_row = pd.read_sql_query(ip_query_list, ip_connection)
            ip_new_row = dict(ip_new_row)
            # print('row:')
            # print(ip_new_row)
            # print(type(ip_new_row))
            ip_new_row = pd.DataFrame(ip_new_row)
            # print(ip_new_row)
            df_pf_dufs = df_pf_dufs.append(ip_new_row, ignore_index=True)
            
        except:
            print(f"Couldn't get data from the Oracle database for {i}")
            list_pf_could_not_get_data = list_pf_could_not_get_data.append(i)

    return df_ip_dufs, df_pf_dufs, list_pf_could_not_get_data

###
# Beatriz Lazaro Pinto (October 2024)
#
# This script produces two tables>
# - InterPro DUFs that are not DUFs in Pfam
# - Pfam DUFs that are not DUFs in InterPro
# with all the information that can help quickly see
# if further curation is worth it
#
# The config_deDUFconsistency.ini file must be in the folder where this script is run
###

if __name__ == "__main__":

    # get information from config file
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f'Cannot open {args.config}: no such file or directory')

    config = ConfigParser()
    config.read(args.config)

    # get pfam info and create dfs

    pf_user = str(config['pfam']['user'])
    pf_pw = str(config['pfam']['password'])
    pf_host = str(config['pfam']['host'])
    pf_port = str(config['pfam']['port'])
    pf_db = str(config['pfam']['db'])

    pf_query = "select pfamA_id,pfamA_acc,description,previous_id from pfamA"

    # connect to Pfam database
    pf_connection = mysql.connector.connect(
        host=pf_host,
        user=pf_user,
        password=pf_pw,
        database=pf_db,
        port=pf_port
    )

    # Check if the connection is successful
    if pf_connection.is_connected():
        print("Succesfully connected to the Pfam database")
        pass
    else:
        print("Failed to connect to the Pfam database")
        exit()

    # create dataframe of all pfam entries
    pf_df = pd.read_sql(pf_query, pf_connection)
    pf_df = pf_df.set_index('pfamA_acc')
    print(pf_df.head())
    pf_connection.close()

    # create dataframe of pfam DUFs
    pf_duf_df = get_pf_duf_df(pf_df) 
    print(pf_duf_df.head(2))

    # get ip info and create df of the DUFs in IP with a Pfam contributing signature

    ip_user = str(config['interpro']['user'])
    ip_pw = str(config['interpro']['password'])
    ip_db = str(config['interpro']['db'])
    
    ip_query = """
        select e.entry_ac, e2m.method_ac, e.name as IP_NAME, e.short_name, em.description, em.name as PF_ID
        from interpro.entry e
        join interpro.entry2method e2m on e.entry_ac=e2m.entry_ac
        join interpro.method em on e2m.method_ac=em.method_ac
        where (e.short_name LIKE '%DUF%'
        or e.name LIKE '%UPF%'
        or e.name LIKE '%UCP%')
        and e.checked='Y'
        and em.dbcode='H'
    """

    try:
        ip_connection = oracledb.connect(
            user=ip_user,
            password=ip_pw,
            dsn=ip_db)
        print("Successfully connected to the InterPro Database.")

        ip_df = pd.read_sql_query(ip_query, ip_connection)
        print(ip_df.head(2))

        try:
            list_mismatches = check_mismatching(ip_df, pf_duf_df)
            print (list_mismatches[:3])

            try:
                df_ip_dufs, df_pf_dufs, list_pf_could_not_get_data = create_df(list_mismatches, ip_df, pf_df) #this function does SQL queries
                
                df_ip_dufs.to_csv('df_ip_dufs_to_check.csv')
                df_pf_dufs.to_csv('df_pf_dufs_to_check.csv')
                
                print('df_ip_dufs created.')
                print('df_pf_dufs created.')
                if list_pf_could_not_get_data:
                    print(list_pf_could_not_get_data)

            except:
                print("Couldn't create the dataframes.")

        except:
            print("Couldn't get the list of mismatches.")
        
        ip_connection.close()

    except:
        print("Couldn't connect to the InterPro database.")
    