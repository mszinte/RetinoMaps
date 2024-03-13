"""
-----------------------------------------------------------------------------------------
170k_gausgridfit_averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
Average all all the subject of the studie on the 170k space.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
2. run python command
>> python 170gaussgridfit_averaging.py [main directory] [project name] [group] 
    [server project]
-----------------------------------------------------------------------------------------
Exemple:
python 170k_gausgridfit_averaging.py /scratch/mszinte/data RetinoMaps 327 
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import os
import sys
import json
import ipdb
import numpy as np
import nibabel as nb
deb = ipdb.set_trace

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import load_surface , make_surface_image

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects_all']


prf_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/prf'.format(main_dir, project_dir)
avg_170k_dir = '{}/prf_derivatives'.format(prf_dir)
os.makedirs(avg_170k_dir, exist_ok=True)



avg_170k_fn = 'sub-170k_task-pRF_fmriprep_dct_avg_prf-deriv_gauss_gridfit.dtseries.nii'
for n_subject, subject in enumerate(subjects) :
    print('adding {} to averaging'.format(subject))
    
    deriv_dir = '{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives'.format(main_dir, project_dir, subject)

    deriv_fn = '{}_task-pRF_fmriprep_dct_avg_prf-deriv_gauss_gridfit.dtseries.nii'.format(subject)

    img, data = load_surface(fn='{}/{}'.format(deriv_dir, deriv_fn))

    
    if n_subject == 0:
        data_avg = np.copy(data)
    else:
        data_avg = np.nanmean(np.array([data_avg, data]), axis=0)

    
#  export results 
maps_names = ['rsq', 'ecc', 'polar_real', 'polar_imag', 'size',
              'amplitude','baseline', 'x','y','hrf_1','hrf_2']

print('saving {}/{}'.format(avg_170k_dir, avg_170k_fn))
avg_img = make_surface_image(data=data_avg, source_img=img, maps_names=maps_names)
nb.save(avg_img,'{}/{}'.format(avg_170k_dir, avg_170k_fn))

# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))