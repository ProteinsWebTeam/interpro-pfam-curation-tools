from datetime import datetime
import os
import subprocess


def prepare_iteration(folder_path):

    '''Functions in a folder (folder_path points to this folder) that contains subfolders with all files related to
    building an HMM model and check overlaps'''

    for file in os.listdir(folder_path):
        pfamout = file + '//PFAMOUT'

        if os.path.isfile(pfamout):
            try:
                members = subprocess.check_output('cd ' + file + '; wc -l ALIGN', shell=True)
                members =str(members)[2:-8]
                with open('smORF_main.log', 'a') as f:
                    current_datetime = str(datetime.now())
                    f.write(f'{current_datetime} For protein sequence: {file} \nMembers: {members} \n')
            except:
                with open('smORF_main.log', 'a') as f:
                    f.write(f'Skipping {file} \n')
                continue
            
            try:
                overlaps = subprocess.check_output('cd ' + file + '; wc -l overlap', shell=True)
                overlaps = str(overlaps)[2:-10]
                with open('smORF_main.log', 'a') as f:
                    f.write('Overlaps: {overlaps} \n')

                if int(members) > 0:    
                    percentage_overlaps = int(overlaps)/int(members)
                    
                    if percentage_overlaps < 0.3:
                        os.system('cd ' + file + '; mkdir iteration; cd iteration; \
                        cp ../ALIGN ./SEED; cp ../DESC .; cp ../SEED ./oldSEED')

                    if percentage_overlaps > 0.9:
                        os.system('mv ' + file + ' overlaps_' + file)

                    with open('smORF_main.log', 'a') as f:
                        f.write(f'Percentage overlaps: {percentage_overlaps} \n\n')
                
                if int(members) == 0:
                    os.system('mv ' + file + ' nomembers_' + file)  
                    with open('smORF_main.log', 'a') as f:
                        f.write(f'Could not calculate the overlap percentage.\n\n')

            except:
                pass

    with open('smORF_main.log', 'a') as f:
        f.write(f'Iterations are ready.\n\n')


                    

    