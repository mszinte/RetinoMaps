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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/glm/
2. run python command
python preproc_end.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python glm_fit.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import nibabel as nb

import numpy as np


import json
import glob
import ipdb


import sys
import os
import importlib
import datetime

# nilearn import

from nilearn.glm import fdr_threshold

from nilearn.glm.first_level import make_first_level_design_matrix,run_glm
from nilearn.glm.contrasts import compute_contrast
import scipy.stats as stats


import warnings
warnings.filterwarnings("ignore")


sys.path.append("{}/../../utils".format(os.getcwd()))
from glm_utils import eventsMatrix
from gifti_utils import make_gifti_image


# Get inputs

start_time = datetime.datetime.now()
deb = ipdb.set_trace

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
TR = analysis_info['TR']
TRs = analysis_info['TRs']
glm_alpha = analysis_info['glm_alpha']
high_pass_type = analysis_info['high_pass_type'] 

session = 'ses-02'

hemis = ['L','R']


tasks = ["pMF","PurLoc","SacLoc","PurVELoc","SacVELoc"]
    

for task in tasks : 
    print(task)
    for hemi in hemis : 
            


        
        # 
        event_dir = '{}/{}/{}/{}/func/'.format(main_dir,project_dir,subject,session)
        event_file = glob.glob("{}/{}_{}_task-{}_run-*_events.tsv".format(event_dir,subject,session,task))
        avg_dir = '{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/fsnative'.format(main_dir,project_dir,subject)
        
        glm_dir = '/{}/{}/derivatives/pp_data/{}/glm'.format(main_dir,project_dir,subject)
        os.makedirs(glm_dir, exist_ok=True)
        
        
        out_pred_name = '{}/{}_hemi-{}_{}_glm_pred.func.gii'.format(glm_dir,subject,hemi,task)
        out_fit_name = '{}/{}_hemi-{}_{}_glm_fit.func.gii'.format(glm_dir,subject,hemi,task)
        
        img_avg_bold_hemi = nb.load('{}/{}_task-{}_hemi-{}_fmriprep_dct_avg_bold.func.gii'.format(avg_dir,subject,task,hemi))
        img_avg_bold_hemi_header = img_avg_bold_hemi.header
        img_avg_bold_hemi_meta = img_avg_bold_hemi.meta
        
        print('data load')
        data_avg_bold_hemi = [x.data for x in img_avg_bold_hemi.darrays]
        data_avg_bold_hemi = np.vstack(data_avg_bold_hemi) 
        
        print(TR)
        events = eventsMatrix(design_file=event_file[1], task=task, tr=TR)
        
        frame_times = np.arange(TRs) * TR
        design_matrix = make_first_level_design_matrix(frame_times,
                                                   events=events,
                                                   hrf_model='spm',
                                                   drift_model=None
                                                   )
        
        if task == 'SacLoc':
            cond1_label, cond2_label = ['Sac'], ['Fix']
            comp_num = 1
        elif task == 'PurLoc':
            cond1_label, cond2_label = ['Pur'], ['Fix']
            comp_num = 1
        elif task == 'SacVELoc':
            cond1_label, cond2_label = ['SacExo','SacExo','SacEndo'], ['SacEndo','Fix','Fix']
            comp_num = 3
        elif task == 'PurVELoc':
            cond1_label, cond2_label = ['PurExo','PurExo','PurEndo'], ['PurEndo','Fix','Fix']
            comp_num = 3
        elif task == 'pMF':
            cond1_label, cond2_label = ['PurSac'], ['Fix']
            comp_num = 1
            
    
        labels_hemi, estimates_hemi = run_glm(data_avg_bold_hemi, design_matrix.values,noise_model="ols")
        
        print('glm done')
        
        pred_hemi = np.zeros(data_avg_bold_hemi.shape)
        rsquare_hemi = np.zeros_like(labels_hemi)
        for label_ in estimates_hemi:
            label_mask = labels_hemi == label_
            reg = estimates_hemi[label_]
            pred_hemi[:,label_mask] = reg.predicted
            rsquare_hemi[label_mask] = reg.r_square
        
        for contrast_num, contrast in enumerate(zip(cond1_label,cond2_label)):
            print(contrast_num)
            
            contrast_values = (design_matrix.columns == contrast[0]) * 1.0 -(design_matrix.columns == contrast[1])
            eff = compute_contrast(labels_hemi, estimates_hemi, contrast_values,contrast_type='t')
        
        
            z_map = eff.z_score()
            z_p_map = 2*(1 - stats.norm.cdf(abs(z_map)))
            fdr_th = fdr_threshold(z_map, glm_alpha)
            fdr = z_map
            fdr *= (z_map > fdr_th)
            fdr_p_map = 2*(1 - stats.norm.cdf(abs(fdr)))
            
            
            
            if contrast_num:
                fit = np.vstack((fit,z_map,z_p_map,fdr,fdr_p_map))
            else:                 
                fit = np.vstack((z_map,z_p_map,fdr,fdr_p_map,rsquare_hemi))
            
        print(fit.shape)        
    
        
    
        
        # export fit param      
        # fit_img_hemi = nb.gifti.GiftiImage(header=img_avg_bold_hemi_header, meta=img_avg_bold_hemi_meta)
    
        # for i in range(fit.shape[0]):
        #     data = fit[i,:]
        #     darray = nb.gifti.GiftiDataArray(data,datatype = 'NIFTI_TYPE_FLOAT32')
        #     fit_img_hemi.add_gifti_data_array(darray)
        # nb.save(fit_img_hemi, out_fit_name)
        # print('export fit done')
            
        # export prediction 
        
        
        fit_img_hemi = make_gifti_image(img_avg_bold_hemi,fit)
        nb.save(fit_img_hemi, out_fit_name)
        
        
        pred_img_hemi = make_gifti_image(img_avg_bold_hemi,pred_hemi)
        nb.save(pred_img_hemi, out_pred_name)
        
        
        maps_names = ['z_map','z_p_map','fdr','fdr_p_map','rsquare_map']
        
        for map_num, mape_name in enumerate(maps_names):
            os.system('wb_command -set-map-names {} -map {} {}'.format(out_fit_name,map_num+1,mape_name))
            

        
    
        
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))

        
        
    
            