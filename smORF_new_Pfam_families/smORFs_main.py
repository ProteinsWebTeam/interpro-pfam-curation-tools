import requests
from configparser import ConfigParser
import argparse
import os
import subprocess
import time
import smORFs_DESC
import smORFs_iteration
from datetime import datetime


# function to get a list of PMIDs of the required papers from the Europepmc API

def get_pmid_list(epmc_query):

    '''Gets a list of PMIDs of the required papers from the EuropePMC API.'''

    epmc_api_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='
    f = '&format=json'

    try:
        response = requests.get(f'{epmc_api_url} {epmc_query} {f}')
        epmc_js = response.json()
        pmid_list = [paper['id'] for paper in epmc_js['resultList']['result']]
        return pmid_list

    except:
        print(f'The EuropePMC API give unexpected results for the query: {epmc_query}')
        with open('smORF_main.log', 'a') as f:
            current_datetime = str(datetime.now())
            f.write(f'{current_datetime} The EuropePMC API give unexpected results for the query: {epmc_query} \n')



# function to get dictionary of lists of protein IDs matching the characteristics for each pmid

def get_dict_pmid_prot(pmid_list, up_query_1, up_query_2):
    
    '''Gets a dictionary of lists of protein IDs matching the characteristics for each PMID.'''

    smORF_dict = {}
    up_api = 'https://rest.uniprot.org/uniprotkb/search?query='
    f = '&format=json'

    try:
        for pmid in pmid_list:
            list_prots = []
            try:
                response = requests.get(up_api + up_query_1 + pmid + up_query_2 + f)
                up_js = response.json()
                list_prots = [protein['primaryAccession'] for protein in up_js['results']]
                smORF_dict[pmid] = list_prots
            except:
                continue

        return smORF_dict

    except:
        print('The Uniprot API give unexpected results')
        with open('smORF_main.log', 'a') as f:
            current_datetime = str(datetime.now())
            f.write(f'{current_datetime} The Uniprot API give unexpected results \n')


# function to build initial folders of potential families with SEEDs
def create_folders_potential_fam(smORF_dict, pfamseq):

    '''Builds folders containing SEEDs for potential new families'''

    for pmid, list_prots in smORF_dict.items():
        with open('smORF_main.log', 'a') as f:
            current_datetime = str(datetime.now())
            f.write(f'{current_datetime} PMID: {pmid} \n')

        for prot in list_prots:
            with open('smORF_main.log', 'a') as f:
                current_datetime = str(datetime.now())
                f.write(f'{current_datetime} Protein sequence: {prot} \n')

            if not os.path.isdir(prot):
                # create folders and initial SEED 
                os.system('mkdir ' + prot + '; cd ' + prot + ';touch ' + pmid + \
                    '; pfetch ' + prot + ' > fasta_SEED && belvu -o Mul fasta_SEED > SEED')
                
                # check if the SEED sequence is in pfamseq
                try:
                    seed = subprocess.check_output(['more ', prot + '/SEED'])
                    seed = str(seed)[2:].split('/')[0]

                    pfseq = str(subprocess.check_output(['esl-sfetch ', pfamseq , seed]))

                    with open('smORF_main.log', 'a') as f:
                        current_datetime = str(datetime.now())
                        f.write(f'{current_datetime} Checking if {seed} is in pfamseq: {pfseq} \n')
                        f.write(f'seed: {seed}')

                except:
                    with open('smORF_main.log', 'a') as f:
                        current_datetime = str(datetime.now())
                        f.write(f'{current_datetime} Could not check if {prot} is in pfamseq \n')

            else:
                with open('smORF_main.log', 'a') as f:
                    f.write(f'Folder was already created for {prot} \n')
            
    with open('smORF_main.log', 'a') as f:
        current_datetime = str(datetime.now())
        f.write(f'{current_datetime} Folders created. \n \n')    
            
def procces_folder(folder_path, author):
            
    '''Functions in a folder (folder_path points to this folder) that contains subfolders with SEEDs 
    ready to be used for the creation of HMMs.'''

    for file in os.listdir(folder_path):
        
        if os.path.isdir(str(file)) or not os.path.isfile(str(file) + '\PFAMOUT'): # check the job was not run already
            try:
                with open('smORF_main.log', 'a') as f:
                    current_datetime = str(datetime.now())
                    f.write(f'\n \n{current_datetime} {file} \n')

                # build model
                cmd = ['pfbuild', '-withpfmake']
                pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=file)
                std_out, std_err = pipes.communicate()

                with open('smORF_main.log', 'a') as f:
                    f.write(f'std_out: {std_out.decode("utf-8")} \n')
                    f.write(f'std_err: {std_err.decode("utf-8")} \n')
                    f.write(f'returncode: {pipes.returncode} \n')

                # if the job was submitted, wait for it to be finished
                if pipes.returncode == 0:

                    if 'LSF_ENVDIR' in os.environ:
                        # with open('smORF_main.log', 'a') as f:
                        #     f.write(f'This is codon. \n')
                        job = std_out.decode("utf-8")[std_out.decode("utf-8").index('<',0)+1:std_out.decode("utf-8").index('>',0)]
                        with open('smORF_main.log', 'a') as f:
                            f.write(f'Job ID: {job} \n')

                        state = subprocess.check_output(['bjobs', '-noheader', job])
                        state = state.decode("utf-8").split(' ')[8]
                        time.sleep(100)

                        with open('smORF_main.log', 'a') as f:
                            f.write(f'The state of the job {job} is: {state} \n')

                        while state != 'DONE':
                            time.sleep(15)
                            state = subprocess.check_output(['bjobs', '-noheader', job])
                            state = state.decode("utf-8").split(' ')[8]
                                    
                        else:
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} The state of the job {job} is: {state} \n')
                    
                    else:
                        # with open('smORF_main.log', 'a') as f:
                        #     f.write('This is Slurm.' + '\n')
                        job = std_out.decode("utf-8").split()[-1]
                        with open('smORF_main.log', 'a') as f:
                            f.write(f'Job ID: {job} \n')

                        time.sleep(5)
                        state = subprocess.check_output(['jobinfo', job])
                        state = state.decode("utf-8").split()[state.decode("utf-8").split().index('State') + 2]
                        time.sleep(95)

                        with open('smORF_main.log', 'a') as f:
                            f.write( f'The state of the job {job} is: {state} \n')

                        while state != 'COMPLETED':                            
                            time.sleep(15)
                            state = subprocess.check_output(['jobinfo', job])
                            state = state.decode("utf-8").split()[state.decode("utf-8").split().index('State') + 2]
                                    
                        else:
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} The state of the job {job} is: {state} \n')

                    if not os.path.isdir(str(file) + '\overlaps'): # Check it has not been processed already

                        # when the job is finished, get SwissProt matches, and add pdb references
                        try:
                            os.system('cd  ' + file + '; swissprot.pl && ~agb/Scripts/add_pdb_ref.pl')
                            sp = subprocess.check_output(['wc','-l', 'sp'], cwd=file)
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} Checking SwissProt matches for {file}.' + '\n')
                        except:
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} Could not get SwissProt matches for {file}. \n')
                        
                        # get taxonomy distribution
                        try:
                            os.system('cd  ' + file + '; python /homes/blp/script_tests/taxonomy.py')
                            with open('smORF_main.log', 'a') as f:
                                    current_datetime = str(datetime.now())
                                    f.write(f'{current_datetime} Getting taxonomy distribution for {file} \n')
                        except:
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} Could not get taxonomy distribution for {file} \n')

                        # check overlaps
                        try:
                            os.system('pqc-overlap-rdb ' + file )
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} Checking overlaps for {file} \n')
                                
                        except:
                            with open('smORF_main.log', 'a') as f:
                                current_datetime = str(datetime.now())
                                f.write(f'{current_datetime} Overlaps could not be checked for {file}. Skipping. \n')

                        # finally, amend the DESC file
                        try:
                            smORFs_DESC.amend_DESC(file, author)
                            with open('smORF_main.log', 'a') as f:
                                f.write(f'Modifying DESC for {file}. \n')
                        except:
                            with open('smORF_main.log', 'a') as f:
                                f.write(f'DESC could not be modified for {file}. \n')
                    else:
                        with open('smORF_main.log', 'a') as f:
                            f.write(f'There is already an overlaps file for {file} . Skipping. \n')
                else:
                    with open('smORF_main.log', 'a') as f:
                        f.write(f'The job for {file} was not submited. \n')
            except:
                pass

        else:
            with open('smORF_main.log', 'a') as f:
                f.write(f'Skipping {file}. \n')
      
    with open('smORF_main.log', 'a') as f:
        current_datetime = str(datetime.now())
        f.write(f'{current_datetime} Finished processing folders. \n\n')


###
# This script gets Uniprot IDs of 
# - small protein sequences 
# - that are not in the release of Pfam in the Uniprot website
# - linked to a EuropePMC query specified in the config.ini file
# and creates potential family folders out of it 
###
    
if __name__ == "__main__":

    print("smORF script running...")

    # get information from config file
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f'Cannot open {args.config}: no such file or directory')

    config = ConfigParser()
    config.read(args.config)

    # create log file
    with open('smORF_main.log', 'w') as f:
        current_datetime_0 = datetime.now()
        f.write(f'Starting time: {str(current_datetime_0)} \n \n')

    # get list of PMIDs
    epmc_query = str(config['epmc']['epmc_query'])
    pmid_list = get_pmid_list(epmc_query)

    try:
        with open('smORF_main.log', 'a') as f:
            f.write(f'List of PMIDs: {pmid_list} \n \n')
    except:
        pass

    # get pertinent Uniprot IDs
    up_query_1 = str(config['uniprot']['up_query_1'])
    up_query_2 = str(config['uniprot']['up_query_2'])
    smORF_dict = get_dict_pmid_prot(pmid_list, up_query_1, up_query_2)
    try:
        with open('smORF_main.log', 'a') as f:
            f.write(f'Dictionary of protein sequences: {smORF_dict} \n \n')
    except:
        pass
           
    # create protential new family folders and record if Uniprot IDs are on pfamseq 
    folder_path = str(config['files']['folder_path'])
    author = str(config['files']['author'])
    pfamseq = str(config['files']['pfamseq'])
    create_folders_potential_fam(smORF_dict,pfamseq)
    procces_folder(folder_path, author)
    smORFs_iteration.prepare_iteration(folder_path)

    # record end time and total time
    with open('smORF_main.log', 'a') as f:
        current_datetime_end = datetime.now()
        f.write(f'Final time: {str(current_datetime_end)} \n \n')
        final_time = (current_datetime_end - current_datetime_0).total_seconds()/60
        f.write(f'Total time employed: {final_time} minutes.')
