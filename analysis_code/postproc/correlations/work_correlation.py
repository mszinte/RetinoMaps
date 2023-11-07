"""
-----------------------------------------------------------------------------------------
preproc_end.py
-----------------------------------------------------------------------------------------
Goal of the script:
High-pass filter, z-score, average data and pick anat files
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# Preprocessed and averaged timeseries files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd /home/mszinte/projects/RetinoMaps/analysis_code/preproc/functional/
2. run python command
python preproc_end.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python preproc_end.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
import json
import sys
import os
import glob
import numpy as np
import nibabel as nb
import itertools as it
from scipy import stats
import matplotlib.pyplot as plt


# # Inputs
# main_dir = sys.argv[1]
# project_dir = sys.argv[2]
# subject = sys.argv[3]
# group = sys.argv[4]

# # load settings
# with open('../../settings.json') as f:
#     json_s = f.read()
#     analysis_info = json.loads(json_s)
# TR = analysis_info['TR']
# high_pass_threshold = analysis_info['high_pass_threshold'] 
# high_pass_type = analysis_info['high_pass_type'] 
# sessions = analysis_info['session']


main_dir = '/Users/uriel/disks/meso_shared'
project_dir = 'RetinoMaps'
subject = 'sub-02'
group = '327'
session = 'ses-01'

# load settings
with open('/Users/uriel/disks/meso_H/projects/RetinoMaps/analysis_code/settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
TR = analysis_info['TR']
#tasks = analysis_info['task_names']
high_pass_threshold = analysis_info['high_pass_threshold'] 
high_pass_type = analysis_info['high_pass_type'] 
sessions = analysis_info['session']



task = 'pRF'

pp_data_func_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/fsnative".format(main_dir, project_dir, subject)

preproc_files = glob.glob("{}/*_task-{}_*_hemi-R_space-fsnative_bold_{}.func.gii".format(pp_data_func_dir,task, high_pass_type))
corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/fsnative".format(main_dir, project_dir, subject)
os.makedirs(corr_dir, exist_ok=True)

combinaisons = list(it.combinations(preproc_files, 2))



cor_file = "{}/{}_task-{}_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,high_pass_type)

img = nb.load(preproc_files[0])
data = [x.data for x in img.darrays]
data = np.vstack(data)
header = img.header
meta = img.meta

task_cor = np.zeros((data.shape[1])) 
cor_final = np.zeros((data.shape[1]))



for t, combi in enumerate(combinaisons):
    
    a_img = nb.load(combi[0])

    
    
    a_data = [x.data for x in a_img.darrays]
    a_data = np.vstack(a_data)
    
    b_img = nb.load(combi[1])

    
    
    b_data = [x.data for x in b_img.darrays]
    b_data = np.vstack(b_data)
    
    
    for vertice in range(a_data.shape[1]):
        
        
        corr, _ = stats.pearsonr(a_data[:,vertice], b_data[:,vertice])

        task_cor[vertice] = corr 
        
        
    cor_final += task_cor/len(combinaisons)
    
    corr_img = nb.gifti.GiftiImage(header = header, meta= meta)
    
    for dta_corr in cor_final:
        darray = nb.gifti.GiftiDataArray(dta_corr)
        corr_img.add_gifti_data_array(darray)
        
    path_des = '/Users/uriel/Desktop/essai.func.gii'
    nb.save(corr_img, path_des)

# pp_data_func_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/HCP_170k".format(main_dir, project_dir, subject)

# preproc_files = glob.glob("{}/*_task-{}_*_space-fsLR_den-170k_bold_{}.dtseries.nii".format(pp_data_func_dir,task, high_pass_type))
# corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/HCP_170k".format(main_dir, project_dir, subject)
# os.makedirs(corr_dir, exist_ok=True)

# combinaisons = list(it.combinations(preproc_files, 2))



# cor_file = "{}/{}_task-{}_fmriprep_{}_correlations_bold.dtseries.nii".format(corr_dir, subject, task,high_pass_type)

# img = nb.load(preproc_files[0])
# task_cor = np.zeros((img.shape[1])) 
# cor_final = np.zeros((img.shape[1]))



# for t, combi in enumerate(combinaisons):
    
#     a = nb.load(combi[0])
#     data_a = a.get_fdata()
    
#     b = nb.load(combi[1])
#     data_b = b.get_fdata()
    
    
#     for vertice in range(data_a.shape[1]):
#         corr, _ = stats.pearsonr(data_a[:,vertice], data_b[:,vertice])

#         task_cor[vertice] = corr 
        
        
#     cor_final += task_cor/len(combinaisons)
  
    
  

  

    
#     cor_img = nb.cifti2.cifti2.Cifti2Image(cor_final)

#     path_des = '/Users/uriel/Desktop'
#     nb.save(cor_img,path_des )
    
# df = np.reshape(cor_final,(1,cor_final.shape[0]))       
        
# cor = stats.pearsonr(data_a,data_b)
# corr_matrix = np.corrcoef(data_a,data_b)

    


