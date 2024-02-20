"""
-----------------------------------------------------------------------------------------
make_sub_all_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
make a sub-all tsv
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# sub-all tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
2. run python command
python make_sub_all_tsv.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python make_sub_all_tsv.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import ipdb
deb = ipdb.set_trace

import pandas as pd 
import sys
import json
import os

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']


# load subjects tsv and concatenate them
data_all = pd.DataFrame()
for subject in subjects :
    print('adding {}'.format(subject))

    tsv_dir ='{}/{}/derivatives/pp_data/{}/fsnative/prf/tsv'.format(main_dir, 
                                                                    project_dir, 
                                                                    subject)
    data = pd.read_table('{}/{}_task-prf_loo.tsv'.format(tsv_dir,subject))
    
    data_all = pd.concat([data_all, data], ignore_index=True)
    
data_all = data_all.rename(columns={'subject': 'sub-origine'})
data_all['subject'] = ['sub-all'] * len(data_all)

# export tsv 
tsv_all_dir = '{}/{}/derivatives/pp_data/sub-all/fsnative/prf/tsv'.format(main_dir, 
                                                                          project_dir)
os.makedirs(tsv_all_dir, exist_ok=True)
    
data_all.to_csv('{}/sub-all_task-prf_loo.tsv'.format(tsv_all_dir), sep="\t", na_rep='NaN',index=False)


# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))









