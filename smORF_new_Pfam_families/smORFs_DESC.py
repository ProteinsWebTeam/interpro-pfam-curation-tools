import os
import subprocess
import uniprotinfo

# function to modify DESC file

def amend_DESC(file,author):

    '''Functions in a folder (folder_path points to this folder) that contains subfolders with a DESC file. 
    Modifies the DESC file, adding an author, adds a ref (PMID stored as file name) and adds a SEED origin (PMID)'''

    try:
        filename = os.fsdecode(file)
        # add pmid
        pmid = subprocess.check_output('cd ' + filename + '; ls -t|tail -n 1', shell=True)
        pmid = str(pmid)[2:-3]
        print(pmid)
        os.system('cd ' + filename + '; add_ref.pl ' + pmid)
        DESC = filename + '//DESC'
        # with open('smORF_main.log', 'a') as f:
        #     f.write('Modifying ' + DESC + '\n')

        information = uniprotinfo.get_up_accession_info(str(filename))
        with open(DESC, 'r+') as f:
            lines = f.readlines()
            new_line = 0

            # add SE line, AU line and CC lines
            for i, line in enumerate(lines):
                if line.startswith('GA'):
                    new_line = i
                    break
            f.seek(0)
            f.truncate()
            if len(information['id']) > 0:
                f.write('ID   ' + str(information['id'])[2:-3] + '\n')
            else:
                f.write('ID   x' + '\n')
            if len(information['description']) > 0:
                f.write('DE   ' + str(information['description'])[2:-3] + '\n')
            else:
                f.write('DE   x' + '\n')

            # add author
            f.write('AU   ' + author + '\n')
            # add SEED origin
            f.write('SE   [1]\n')
            
            for line in range(new_line,len(lines)):
                f.write(lines[line])

            # add CC lines
            f.write('CC   This protein family includes\nCC   and similar sequences from\nCC   [1].\n')

            # f.write('CC   ' + str(information['cc'])[2:-3])

    except:
        pass
