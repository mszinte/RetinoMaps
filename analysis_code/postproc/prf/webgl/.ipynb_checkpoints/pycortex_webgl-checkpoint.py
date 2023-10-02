"""
-----------------------------------------------------------------------------------------
pycortex_webgl.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create webgl pycortex across tasks
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: recache db (1 = True, 0 = False)
sys.argv[5]: send to invibe server (1 = True, 0 = False)
--------------------------------------------------------------------------------------- --
Output(s):
pycortex webgl per subject
-----------------------------------------------------------------------------------------
To run: 
0. ON MESO SERVER
1. cd to function
>> cd ~/projects/amblyo_prf/analysis_code/postproc/prf/webgl
2. run python command
>> python pycortex_webgl.py [main directory] [project name] [subject num] [recache] [invibe]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_webgl.py /scratch/mszinte/data/ amblyo_prf sub-01 1 1 
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import cortex
import sys
import importlib
import json
import numpy as np
import matplotlib.pyplot as plt
import ipdb
import h5py
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
deb = ipdb.set_trace

# Get inputs
# ----------
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
recache = bool(int(sys.argv[4]))
webapp = bool(int(sys.argv[5]))

# Define analysis parameters
# --------------------------
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
task = analysis_info["task"]

# create folder
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
datasets_dir = '{}/{}/derivatives/pp_data/{}/prf/pycortex/datasets'.format(main_dir, project_dir, subject)
webgl_dir = '{}/{}/derivatives/webgl/{}/'.format(main_dir, project_dir, subject)
os.makedirs(webgl_dir, exist_ok=True)

# Set pycortex db and colormaps
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# Load datasets
avg_dataset_fn = "{}/{}_task-{}_avg.hdf".format(datasets_dir, subject, task)
avg_dataset = cortex.load(avg_dataset_fn)
loo_avg_dataset_fn = "{}/{}_task-{}_loo_avg.hdf".format(datasets_dir, subject, task)
loo_avg_dataset = cortex.load(loo_avg_dataset_fn)
new_dataset = cortex.Dataset(avg_dataset=avg_dataset, loo_avg_dataset=loo_avg_dataset)
cortex.webgl.make_static(outpath=webgl_dir, data=new_dataset, recache=recache)

# # Send to webapp
# # --------------
# if webapp == True:
#     webapp_dir = '{}{}_{}/'.format(analysis_info['webapp_dir'], subject, preproc)
#     os.system('rsync -avuz --progress {local_dir} {webapp_dir}'.format(local_dir=webgl_dir, webapp_dir=webapp_dir))
#     print('go to : https://invibe.nohost.me/amblyo_prf/{}_{}'.format(subject, preproc))