"""
-----------------------------------------------------------------------------------------
pycortex_import.py
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
python pycortex_import.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Executions:
python pycortex_import.py /scratch/mszinte/data amblyo_prf sub-01 327
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
import json
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

# define directories and get fns
fmriprep_dir = "{}/{}/derivatives/fmriprep".format(main_dir, project_dir)
#jobs_dir = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
#fs_license = "{}/{}/code/freesurfer/license.txt".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
temp_dir = "{}/{}/derivatives/temp_data/{}_rand_ds/".format(main_dir, project_dir, subject)
file_list = sorted(glob.glob("{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/*{}*.nii.gz".format(main_dir, project_dir, subject, task)))
#os.makedirs(jobs_dir, exist_ok=True)

# # Define freesurfer cmd
# freesurfer_cmd = """\
# export FREESURFER_HOME={}/{}/code/freesurfer
# export SUBJECTS_DIR={}\n\
# export FS_LICENSE={}\n\
# source $FREESURFER_HOME/SetUpFreeSurfer.sh\n""".format(main_dir, project_dir, fs_dir, fs_license)


# subprocess.run(["bash", "-c", freesurfer_cmd])





# # define freesurfer command
# freesurfer_cmd = """\
# export FREESURFER_HOME={}/{}/code/freesurfer
# export SUBJECTS_DIR={}\n\
# export FS_LICENSE={}\n\
# source $FREESURFER_HOME/SetUpFreeSurfer.sh\n""".format(main_dir, project_dir, fs_dir, fs_license)

# # define permission cmd
# chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir)
# chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir, group=group)


# # create sh folder and file
# import_freesurfer = "{}/import_freesurfer.sh".format(jobs_dir)

# of = open(import_freesurfer, 'w')
# of.write("{}{}{}".format(freesurfer_cmd,chmod_cmd,chgrp_cmd))
# of.close()


# subprocess.run("import_freesurfer.sh")

# set pycortex db and colormaps
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# add participant to pycortex db
print('import subject in pycortex')
cortex.freesurfer.import_subj(subject, subject, fs_dir, 'smoothwm')

# add participant flat maps
print('import subject flatmaps')
try: cortex.freesurfer.import_flat(fs_subject=subject, cx_subject=subject, 
                                  freesurfer_subject_dir=fs_dir, patch='full', auto_overwrite=True)
except: pass

# add transform to pycortex db
print('Add transform: xfm_name: {}'.format(xfm_name))
transform = cortex.xfm.Transform(np.identity(4), file_list[0])
transform.save(subject, xfm_name, 'magnet')

# add masks to pycortex transform
print('add mask: xfm_name: {}'.format(xfm_name))
xfm_masks = analysis_info['xfm_masks']
ref = nb.load(file_list[0])
for xfm_mask in xfm_masks:
    mask = cortex.get_cortical_mask(subject=subject, xfmname=xfm_name, type=xfm_mask)
    mask_img = nb.Nifti1Image(dataobj=mask.transpose((2,1,0)), affine=ref.affine, header=ref.header)
    mask_file = "{}/db/{}/transforms/{}/mask_{}.nii.gz".format(cortex_dir, subject, xfm_name, xfm_mask)
    mask_img.to_filename(mask_file)

# create participant pycortex overlays
print('create subject pycortex overlays to check')
voxel_vol = cortex.Volume(np.random.randn(ref.shape[2], ref.shape[1], ref.shape[0]), subject = subject, xfmname = xfm_name)
ds = cortex.Dataset(rand=voxel_vol)
cortex.webgl.make_static(outpath=temp_dir, data=ds)

# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))