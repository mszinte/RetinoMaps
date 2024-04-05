"""
-----------------------------------------------------------------------------------------
compute_glm_derivatives.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute glm derivatives 
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/glm/
2. run python command
>> python compute_glm_derivatives.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
python compute_glm_derivatives.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import sys
import glob
import ipdb
import json
import numpy as np
import nibabel as nb
import scipy.stats as stats

# nilearn import
from nilearn.glm import fdr_threshold

# Personal import 
sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface

deb = ipdb.set_trace


# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
tasks = analysis_info['task_glm']
glm_alpha = analysis_info['glm_alpha']


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

z_map_idx, fdr_p_map_idx, r_map_idx = 0, 3, 4


pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
for format_, extension in zip(formats, extensions):
    # Define directories
    glm_deriv_dir = "{}/{}/{}/glm/glm_derivatives".format(pp_dir, subject, format_)
    os.makedirs(glm_deriv_dir, exist_ok=True)
        
    
# find all the glm files 
glm_fit_fns = []
for format_, extension in zip(formats, extensions):
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/glm_fit/*fit*.{}".format(main_dir, 
                                                                                    project_dir, 
                                                                                    subject, 
                                                                                    format_, 
                                                                                    extension))
    glm_fit_fns.extend(list_)
            
# split filtered files  depending of their nature
glm_fsnative_hemi_L, glm_fsnative_hemi_R, glm_170k = [], [], []
for subtype in glm_fit_fns:
    if "hemi-L" in subtype:
        glm_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        glm_fsnative_hemi_R.append(subtype)
    elif "170k" in subtype:
        glm_170k.append(subtype)
        
glm_files_list = [glm_fsnative_hemi_L, 
                  glm_fsnative_hemi_R, 
                  glm_170k]


# Averaging
for glm_files in glm_files_list:
    
    img, data = load_surface(fn=glm_files[0])
    final_map = np.zeros((8,data[fdr_p_map_idx].shape[0]))
    # final_map = np.full((1, data[fdr_p_map_idx].shape[0]), np.nan)
    if glm_files[0].find('hemi-L') != -1: hemi = 'hemi-L'
    elif glm_files[0].find('hemi-R') != -1: hemi = 'hemi-R'
    else: hemi = None
    
    for task in tasks:
        # defind output files names 
        glm_files_task = [file for file in glm_files if task in file]


        # Averaging computation
        glm_img, glm_data = load_surface(fn=glm_files_task[0])
        z_map_avg = np.zeros(glm_data[z_map_idx].shape)
        r_map_avg = np.zeros(glm_data[r_map_idx].shape)
        
        
        for glm_file in glm_files_task:
            glm_img, glm_data = load_surface(fn=glm_file)
            z_map_avg += glm_data[z_map_idx,:]/len(glm_files_task)
            r_map_avg += glm_data[r_map_idx,:]/len(glm_files_task)


        z_p_map = 2*(1 - stats.norm.cdf(abs(z_map_avg)))
        fdr_th = fdr_threshold(z_map_avg, glm_alpha)
        fdr = z_map_avg
        fdr *= (z_map_avg > fdr_th)
        fdr_p_map = 2*(1 - stats.norm.cdf(abs(fdr)))
        
        fit_avg = np.vstack((z_map_avg, z_p_map, fdr, fdr_p_map, r_map_avg))


        # export averaged data
        if hemi:
            avg_fn = "{}/{}/derivatives/pp_data/{}/fsnative/glm/glm_derivatives/{}_task-{}_{}_space-fsnative_dct_glm-fit_avg.func.gii".format(
                main_dir, project_dir, subject, subject, task, hemi)
            
        else:
            avg_fn = "{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/{}_task-{}_space-fsLR_den-170k_dct_glm-fit_avg.dtseries.nii".format(
                main_dir, project_dir, subject, subject, task)

        print('avg save: {}'.format(avg_fn))
        maps_names = ['z_map','z_p_map','fdr','fdr_p_map', 'r_map']
        avg_img = make_surface_image(data=fit_avg, source_img=glm_img, maps_names=maps_names)
        nb.save(avg_img, avg_fn)
        
        # make a final map with both the tasks 

        if task == 'PurLoc':
            task_idx = 1

        elif task == 'SacLoc':
            task_idx = 2
            
        elif task == 'pRF':
            task_idx = 4
        
        for vert, fdr_value in enumerate(fdr_p_map):
            if fdr_value < glm_alpha:
                final_map[0,vert] += task_idx
        

            

    # export final map
    if hemi:
        final_fn = "{}/{}/derivatives/pp_data/{}/fsnative/glm/glm_derivatives/{}_task-eyes-mvt_{}_space-fsnative_dct_glm-significant_map.func.gii".format(
            main_dir, project_dir, subject, subject, hemi)
        
    else:
        final_fn = "{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/{}_task-eyes-mvt_space-fsLR_den-170k_dct_glm-significant_map.dtseries.nii".format(
            main_dir, project_dir, subject, subject)
    

    
    #  Make specifique maps 
    for vert, final_value in enumerate(final_map[0,:]):
        if final_value == 1 or final_value == 3:
            final_map[1,vert] = 1
            
        elif final_value == 2 or final_value == 3:
            final_map[2,vert] = 1
            
        elif final_value == 3 :
            final_map[3,vert] = 1
            
        elif final_value == 4 :
            final_map[4,vert] = 1
            
        elif final_value == 5 :
            final_map[5,vert] = 1
            
        elif final_value == 6 :
            final_map[6,vert] = 1
            
        elif final_value == 7 :
            final_map[7,vert] = 1
        
    
    maps_names_2 = ['all','pursuit','saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuite', 'vision_and_saccade', 'vision_and_pursuite_and_pursuit']
    final_img = make_surface_image(data=final_map, source_img=img, maps_names=maps_names_2)
    nb.save(final_img, final_fn)

    

# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))

