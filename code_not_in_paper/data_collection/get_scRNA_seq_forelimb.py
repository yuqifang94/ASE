#!/usr/bin/env python3.8
import requests
import json
import subprocess
import os
all_datasets = {"E10_5":["ENCSR874BOF","ENCSR874BOF"],
                "E14":["ENCSR792LNZ"],
                "E13":["ENCSR156HHE","ENCSR476JWN"],
                "E15":["ENCSR306UHC","ENCSR337QTQ"],
                "E12":["ENCSR745DNK"],
                "E11":["ENCSR642XVO"],
                "E13_5":["ENCSR787QXE"]
}
HEADERS = {'accept': 'application/json'}
outputDir="../downstream/data/mouse_scRNA/fastq_files/"
os.system("mkdir -p "+outputDir)
for age in all_datasets.keys():
    for accession in all_datasets[age]:
        URL = "https://www.encodeproject.org/reference-epigenomes/"+accession+"/?frame=embedded"
        response = requests.get(URL, headers=HEADERS)
        response_json_dict = response.json()
        sampleFiles=response_json_dict[u'files']
        sampleOutputDir="../downstream/data/mouse_scRNA/fastq_files/"+age+"_"+accession+"/"  
        os.system("mkdir -p "+sampleOutputDir) 
        for fileJSON in sampleFiles:
            sampleAc=fileJSON['accession']
            fastqName=fileJSON['submitted_file_name'].split('/')[-1]
            output_filename=sampleOutputDir+fastqName
            url = 'https://www.encodeproject.org'+fileJSON['href']
            print(sampleAc)
            print(url)
            print(output_filename)
            
            cmd_dl=['curl', '--keepalive-time', '2','-RL',url,"-o",output_filename]
            subprocess.check_call(cmd_dl)

