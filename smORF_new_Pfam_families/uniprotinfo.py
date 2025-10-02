import requests
from datetime import datetime

def get_up_accession_info(up_accession):
    information = {}
    information['description'] = []
    information['id'] = []
    information['cc'] = []

    try:
        up_api ='https://www.ebi.ac.uk/proteins/api/proteins/'
        response = requests.get(up_api + up_accession)
        result = response.json()
        try:
            description = result['protein']['recommendedName']['fullName']
            if type(description) == str:
                information['description'].append(description)
        except:
            pass            
        try:
            description2 = result['protein']['recommendedName']['fullName']['value']
            if type(description2) == str:
                information['description'].append(description2)
        except:
            pass
        try:
            short_name = result['protein']['recommendedName']['shortName'][0]['value']
            information['id'].append(short_name)

        except:
            pass

        try:
            gene = result['gene'][0]['name']['value']
            information['id'].append(gene)
        except:
            pass

        try:
            function = result['comments'][0]['text'][0]['value']
        
            information['cc'].append(function)
        except:
            pass
        
    except:
        with open('smORF_main.log', 'a') as f:
            current_datetime = str(datetime.now())
            f.write(f'{current_datetime} Could not get any information for {up_accession} \n')

    return information


