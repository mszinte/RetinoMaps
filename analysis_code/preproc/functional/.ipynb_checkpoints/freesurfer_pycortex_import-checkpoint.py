"""
-----------------------------------------------------------------------------------------
freesurfer_pycortex_import.py
-----------------------------------------------------------------------------------------
Goal of the script:
Import subject in pycortex database
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/stereo_prf/analysis_code/preproc/functional/
2. run python command
python freesurfer_pycortex_import.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Executions:
python freesurfer_pycortex_import.py /scratch/mszinte/data amblyo_prf sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# imports
import os
import sys
import glob
import numpy as np
import ipdb
import platform
import importlib
import cortex
import nibabel as nb
import subprocess 
deb = ipdb.set_trace

# functions import
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file

# get inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
xfm_name = analysis_info['xfm_name']
task = analysis_info['task']


jobs_dir = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
fs_license = "{}/{}/code/freesurfer/license.txt".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
os.makedirs(jobs_dir, exist_ok=True)

# define freesurfer command
freesurfer_cmd = """\
export FREESURFER_HOME={}/{}/code/freesurfer
export SUBJECTS_DIR={}\n\
export FS_LICENSE={}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh\n""".format(main_dir, project_dir, fs_dir, fs_license)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir, group=group)

#define pycortex cmd
py_cortex_cmd = "python pycortex_import.py {} {} {} {}".format(main_dir,project_dir,subject,group)

# create sh folder and file
import_freesurfer = "{}/import_freesurfer.sh".format(jobs_dir)

of = open(import_freesurfer, 'w')
of.write("{}{}{}{}".format(freesurfer_cmd,chmod_cmd,chgrp_cmd,py_cortex_cmd))
of.close()


os.system("{}".format(import_freesurfer))