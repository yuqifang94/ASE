#!/usr/bin/env python2.7
import requests, json
import subprocess
import os

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}
#all_datasets = ['ENCSR099RYD',
#		'ENCSR916GKL',
#		'ENCSR254OTS',
#		'ENCSR377CZO',
#		'ENCSR700FCF',
#		'ENCSR695BII',
#		'ENCSR473ZCM',
#		'ENCSR950OMB',
#		'ENCSR993XPD',
#		'ENCSR735PMJ',
#		'ENCSR587ASA',
#		'ENCSR628KWF',
#		'ENCSR144UTF',
#		'ENCSR450UMG',
#		'ENCSR613BMI',
#		'ENCSR476OEL',
#		'ENCSR651PAZ',
#		'ENCSR217TMK',
#		'ENCSR089FFK',
#		'ENCSR353IFP',
#		'ENCSR842QTB',
all_datasets =[	'ENCSR888JFA',
		'ENCSR007JSP',
		'ENCSR334GBD',
		'ENCSR835OJU',
		'ENCSR835OJU',
		'ENCSR609ZCD',
		'ENCSR324XQF',
		'ENCSR634OPL',
		'ENCSR803ICQ']
all_datasets=['ENCSR258MDR']

#os.system("mkdir -p CpG_met/")

for accession in all_datasets:
    URL = "https://www.encodeproject.org/reference-epigenomes/"+accession+"/?frame=embedded"

	# GET the object
    response = requests.get(URL, headers=HEADERS)

    # Extract the JSON response as a python dict
    response_json_dict = response.json()
    tissue=response_json_dict['biosample_ontology']['term_name']
    stage=response_json_dict['aliases'][0].split(':')[1].split('_')[0]
    for ind_related_dataset in range(len(response_json_dict['files'])):
        related_dataset = response_json_dict["files"][ind_related_dataset]
        # Get information
        if related_dataset['output_type'] != "methylation state at CpG":
           continue
        if related_dataset['file_type'] !="bed bedMethyl":
	   continue
                                 
        # Download the fastq files and label them with biological replication number
        bio_rep = related_dataset['biological_replicates']
	if len(bio_rep)>1:
	    continue
        output_prefix = "_".join([stage,tissue,str(bio_rep[0])])
        #print output_prefix
        # Download
        url = 'https://www.encodeproject.org'+related_dataset['href']
        print(url)
	output_filename = "CpG_met/" + output_prefix  + ".bed.gz"
        subprocess.check_call(['curl',
                                '-RL',
                                 url,
                                 "-o",
                                 output_filename])
