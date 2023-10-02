"""
-----------------------------------------------------------------------------------------
bad_run.py
-----------------------------------------------------------------------------------------
Goal of the script:
Exclude bad run
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
----------------------------------------------------
Output(s):
re right the bad run names with _exluded extention. 
Supposed to have change the 'setting.json'.
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/disks/meso_H/projects/stereo_prf/analysis_code/preproc/functional
2. run python command
python bad_run.py [main directory] [project name] [subject num]
-----------------------------------------------------------------------------------------
Exemple:
python bad_run.py ~/disks/meso_shared amblyo_prf
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""
# General imports
import json
import os
import sys

# Inputs
main_dir = sys.argv[1]
project_name = sys.argv[2]

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
session = analysis_info['session']
subject_exluded = analysis_info['subject_exluded']
run_exluded = analysis_info['run_exluded']
exclusion_nb = len(run_exluded)

# add the _excluded extention to the bad_run 
for t in range(exclusion_nb):
    folder_path = "{main_dir}/{project_name}/derivatives/fmriprep/fmriprep/{subject_exluded}/ses-02/func".format(main_dir=main_dir,project_name=project_name,subject_exluded=subject_exluded[t])  
    prefix = "{subject_exluded}_ses-02_task-prf_{run_exluded}_".format(subject_exluded=subject_exluded[t],run_exluded=run_exluded[t]  

    for filename in os.listdir(folder_path):
        if filename.startswith(prefix):
            new_filename = filename + "_exluded"
            os.rename(os.path.join(folder_path, filename), os.path.join(folder_path, new_filename))




