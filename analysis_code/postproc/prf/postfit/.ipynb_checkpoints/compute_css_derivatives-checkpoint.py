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
python compute_css_derivatives.py /scratch/mszinte/data RetinoMaps sub-09 327
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
import pandas as pd
import numpy as np
import nibabel as nb


sys.path.append("{}/../../../utils".format(os.getcwd()))
from prf_utils import fit2deriv
from surface_utils import make_surface_image , load_surface
from pycortex_utils import get_roi_masks_hemi
deb = ipdb.set_trace


# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Define directories
prf_deriv_dir = "{}/{}/fsnative/prf/prf_derivatives".format(pp_dir, subject)
os.makedirs(prf_deriv_dir, exist_ok=True)

prf_tsv_dir = "{}/{}/fsnative/prf/tsv".format(pp_dir, subject)
os.makedirs(prf_tsv_dir, exist_ok=True)

# Get prf fit filenames
fit_fns= glob.glob("{}/{}/fsnative/prf/fit/*prf-fit_css*".format(pp_dir,subject))


# Compute derivatives
for fit_fn in fit_fns:
    
    deriv_fn = fit_fn.split('/')[-1].replace('prf-fit', 'prf-deriv')


    if os.path.isfile(fit_fn) == False:
        sys.exit('Missing files, analysis stopped : {}'.format(fit_fn))
    else:
        print('Computing derivatives: {}'.format(deriv_fn))
        
        # get arrays
        fit_img, fit_data = load_surface(fit_fn)

        

 
        # compute and save derivatives array
        maps_names = ['prf_rsq', 'prf_ecc', 'polar_real', 'polar_imag', 'prf_size',
                      'amplitude','baseline', 'prf_x','prf_y','hrf_1','hrf_2','prf_n','prf_loo_r2']
        
        deriv_array = fit2deriv(fit_array=fit_data,model='css',is_loo_r2=True)
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



df_rois_brain = pd.DataFrame()
# Averaging
for loo_deriv_fns in loo_deriv_fns_list:

    # Averaging computation
    deriv_img, deriv_data = load_surface(fn=loo_deriv_fns[0])
    loo_deriv_data_avg = np.zeros(deriv_data.shape)
    for loo_deriv_fn in loo_deriv_fns:
        loo_deriv_avg_fn = loo_deriv_fn.split('/')[-1]
        loo_deriv_avg_fn = re.sub(r'avg_loo-\d+_prf-deriv', 'prf-deriv-loo-avg', loo_deriv_avg_fn)
        
        # load data 
        loo_deriv_img, loo_deriv_data = load_surface(fn=loo_deriv_fn)
        
        
        # Averagin
        loo_deriv_data_avg += loo_deriv_data/len(loo_deriv_fns)
    
    
    # export averaged data in surface format 
    loo_deriv_img = make_surface_image(data=loo_deriv_data_avg, source_img=loo_deriv_img, maps_names=maps_names)
    nb.save(loo_deriv_img,'{}/{}'.format(prf_deriv_dir,loo_deriv_avg_fn))
    
    
    # make an final df for tsv exportation
    roi_verts, hemi = get_roi_masks_hemi(fn=loo_deriv_fns[0], 
                                     subject=subject,
                                     rois=rois)
    
    df_rois_hemi = pd.DataFrame()
    for roi in roi_verts.keys():
        # make a dictionaire with data for each rois 
        data_dict = {col: loo_deriv_data_avg[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
        data_dict['rois'] = [roi] * loo_deriv_data_avg[:, roi_verts[roi]].shape[1]
        data_dict['subject'] = [subject] * loo_deriv_data_avg[:, roi_verts[roi]].shape[1]
        data_dict['hemi'] = [hemi] * loo_deriv_data_avg[:, roi_verts[roi]].shape[1]
    
        
        # make the final dataframe for one hemi
        df_roi = pd.DataFrame(data_dict)
        df_rois_hemi = pd.concat([df_rois_hemi, df_roi], ignore_index=True)

    # make the final dataframe for the brain 
    df_rois_brain = pd.concat([df_rois_brain, df_rois_hemi], ignore_index=True)
    
    
    
# export averaged data in tsv 
df_rois_brain.to_csv('{}/{}_task-prf_loo.tsv'.format(prf_tsv_dir,subject), sep="\t", na_rep='NaN',index=False)
          


# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))

