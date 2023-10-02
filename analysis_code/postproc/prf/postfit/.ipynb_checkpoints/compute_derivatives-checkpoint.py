"""
-----------------------------------------------------------------------------------------
compute_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute pRF derivatives
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/amblyo_prf/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_derivatives.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
python compute_derivatives.py /scratch/mszinte/data amblyo_prf sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import glob
import ipdb
import json
import nibabel as nb
import numpy as np
import os
import sys
sys.path.append("{}/../../../utils".format(os.getcwd()))
from prf_utils import fit2deriv
deb = ipdb.set_trace

# Define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
task = analysis_info['task']

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
prf_fit_dir = "{}/{}/prf/fit".format(pp_dir, subject)

# Get timeseries filenames
pp_avg_fns = glob.glob("{}/{}/func/fmriprep_dct_avg/*avg*.nii.gz".format(pp_dir,subject))

# Compute derivatives
for pp_avg_fn in pp_avg_fns:
    
    fit_fn = "{}/{}_prf-fit.nii.gz".format(prf_fit_dir, os.path.basename(pp_avg_fn)[:-7])
    pred_fn = "{}/{}_prf-pred.nii.gz".format(prf_fit_dir, os.path.basename(pp_avg_fn)[:-7])
    deriv_fn = "{}/{}_prf-deriv.nii.gz".format(prf_fit_dir, os.path.basename(pp_avg_fn)[:-7])
    
    if os.path.isfile(fit_fn) == False:
        sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
    else:
        print('Computing derivatives: {}'.format(deriv_fn))
        
        # get arrays
        fit_img = nb.load(fit_fn)
        fit_array = fit_img.get_fdata()
        data_array = nb.load(pp_avg_fn).get_fdata()
        pred_array = nb.load(pred_fn).get_fdata()
        
        # compute and save derivatives array
        deriv_array = fit2deriv(fit_array=fit_array, data_array=data_array, pred_array=pred_array)
        deriv_img = nb.Nifti1Image(dataobj=deriv_array, affine=fit_img.affine, header=fit_img.header)
        deriv_img.to_filename(deriv_fn)

# compute average loo derivatives
loo_deriv_avg_fn = "{}/{}_task-{}_fmriprep_dct_bold_loo_avg_prf-deriv.nii.gz".format(prf_fit_dir, subject, task)
print('Computing derivatives: {}'.format(loo_deriv_avg_fn))

loo_deriv_fns = glob.glob("{}/*loo*{}-deriv.nii.gz".format(prf_fit_dir, task))
loo_deriv_array = np.zeros_like(deriv_array)
for loo_deriv_fn in loo_deriv_fns:
    loo_deriv_array += nb.load(loo_deriv_fn).get_fdata()/len(loo_deriv_fns)
loo_deriv_img = nb.Nifti1Image(dataobj=loo_deriv_array, affine=fit_img.affine, header=fit_img.header)
loo_deriv_img.to_filename(loo_deriv_avg_fn)

# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))