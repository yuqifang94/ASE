#!/usr/bin/env python2.7
import requests, json
import subprocess
import os

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}
all_datasets = ['ENCSR839HKP',
		'ENCSR790SUM',
		'ENCSR296HYS',
		'ENCSR082AVR',
		'ENCSR369AYR',
		'ENCSR781CPK',
		'ENCSR414HVF',
		'ENCSR242GOA']


os.system("mkdir -p CpG_met/")

for tissue_accession in all_datasets:
    URL_ts= "https://www.encodeproject.org/reference-epigenomes/"+tissue_accession+"/?frame=embedded"
    response_ts=requests.get(URL_ts,headers=HEADERS)
    response_json_dict_ts=response_ts.json()
    for stage_idx in range(len(response_json_dict_ts['related_datasets'])):
        accession=response_json_dict_ts['related_datasets'][stage_idx]['accession']
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
