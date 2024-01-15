"""
-----------------------------------------------------------------------------------------
compute_css_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute pRF derivatives from the pRF grid gauss fit
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_derivatives.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import re
import sys
import glob
import ipdb
import json
import numpy as np
import nibabel as nb


sys.path.append("{}/../../../utils".format(os.getcwd()))


from prf_utils import fit2deriv
from surface_utils import make_surface_image , load_surface
deb = ipdb.set_trace


# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Define directories
prf_fit_dir = "{}/{}/fsnative/prf/fit".format(pp_dir, subject)
prf_deriv_dir = "{}/{}/fsnative/prf/prf_derivatives".format(pp_dir, subject)
os.makedirs(prf_deriv_dir, exist_ok=True)

# Get prf fit filenames
fit_fns= glob.glob("{}/{}/fsnative/prf/fit/*prf-fit_css*".format(pp_dir,subject))

# Compute derivatives
for fit_fn in fit_fns:
    
    deriv_fn = fit_fn.split('/')[-1]
    deriv_fn = deriv_fn.replace('prf-fit', 'prf-deriv')


    if os.path.isfile(fit_fn) == False:
        sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
    else:
        print('Computing derivatives: {}'.format(deriv_fn))
        
        # get arrays
        fit_img, fit_data = load_surface(fit_fn)
        

 
        # compute and save derivatives array
        maps_names = ['rsq', 'ecc', 'polar_real', 'polar_imag', 'size',
                      'amplitude','baseline', 'x','y','hrf_1','hrf_2','n']
        
        deriv_array = fit2deriv(fit_array=fit_data,model='css')
        deriv_img = make_surface_image(data=deriv_array, source_img=fit_img, maps_names=maps_names)
        nb.save(deriv_img,'{}/{}'.format(prf_deriv_dir,deriv_fn))


# compute prf derivatives average of loo 

# find all the filtered files 
derives_fns = glob.glob("{}/*loo-*_prf-deriv_css.func.gii".format(prf_deriv_dir))
            


# split filtered files  depending of their nature
deriv_fsnative_hemi_L, deriv_fsnative_hemi_R = [], []
for subtype in derives_fns:
    if "hemi-L" in subtype:
        deriv_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        deriv_fsnative_hemi_R.append(subtype)

        
loo_deriv_fns_list = [deriv_fsnative_hemi_L,
                      deriv_fsnative_hemi_R]



# Averaging
for loo_deriv_fns in loo_deriv_fns_list:

    # defind output files names 

    if loo_deriv_fns[0].find('hemi-L') != -1: hemi = 'hemi-L'
    elif loo_deriv_fns[0].find('hemi-R') != -1: hemi = 'hemi-R'


    # Averaging computation
    deriv_img, deriv_data = load_surface(fn=loo_deriv_fns[0])
    loo_deriv_data_avg = np.zeros(deriv_data.shape)
    for loo_deriv_fn in loo_deriv_fns:
        loo_deriv_avg_fn = loo_deriv_fn.split('/')[-1]
        loo_deriv_avg_fn = re.sub(r'avg_loo-\d+_prf-deriv', 'prf-deriv-loo-avg', loo_deriv_avg_fn)

        
        loo_deriv_img, loo_deriv_data = load_surface(fn=loo_deriv_fn)
        loo_deriv_data_avg += loo_deriv_data/len(loo_deriv_fns)

    # export averaged data
    loo_deriv_img = make_surface_image(data=loo_deriv_data_avg, source_img=loo_deriv_img, maps_names=maps_names)
    nb.save(loo_deriv_img,'{}/{}'.format(prf_deriv_dir,loo_deriv_avg_fn))
          


# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))

