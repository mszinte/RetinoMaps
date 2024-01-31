"""
-----------------------------------------------------------------------------------------
glm_fit_concat.py
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
python glm_fit_concat.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python glm_fit_concat.py /scratch/mszinte/data RetinoMaps sub-07 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug input
import ipdb 
deb = ipdb.set_trace

# General imports
import os
import re
import sys
import json
import glob
import datetime
import numpy as np
import pandas as pd
import nibabel as nb
import scipy.stats as stats

# nilearn import
from nilearn.glm import fdr_threshold
from nilearn.glm.contrasts import compute_contrast
from nilearn.glm.first_level import make_first_level_design_matrix, run_glm

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from glm_utils import eventsMatrix, extract_predictions_r2
from surface_utils import load_surface, make_surface_image


# Start counting the elapsed time for code execution
start_time = datetime.datetime.now()


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
tasks = analysis_info['task_glm']
func_session = analysis_info['func_session'][0]
formats = analysis_info['formats']
extensions = analysis_info['extensions']
confounds_list = analysis_info['glm_confounds']






for format_, extension in zip(formats, extensions):
    # make folders
    glm_dir = '/{}/{}/derivatives/pp_data/{}/{}/glm/glm_fit'.format(main_dir, 
                                                                    project_dir, 
                                                                    subject, 
                                                                    format_)
    os.makedirs(glm_dir, exist_ok=True)
    
preproc_fns = []
for format_, extension in zip(formats, extensions):
    list_ = glob.glob('{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct/*dct*.{}'.format(main_dir, 
                                                                                                        project_dir, 
                                                                                                        subject, 
                                                                                                        format_,  
                                                                                                       extension))
    preproc_fns.extend(list_)



# split filtered files  depending of their nature
preproc_fsnative_hemi_L, preproc_fsnative_hemi_R, preproc_170k = [], [], []
for subtype in preproc_fns:
    if "hemi-L" in subtype:
        preproc_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        preproc_fsnative_hemi_R.append(subtype)
    elif "170k" in subtype:
        preproc_170k.append(subtype)
            
preproc_files_list = [preproc_fsnative_hemi_L, 
                      preproc_fsnative_hemi_R, 
                      preproc_170k]

for preproc_fn_runs in preproc_files_list :
    if preproc_fn_runs[0].find('hemi-L') != -1: hemi = 'hemi-L'
    elif preproc_fn_runs[0].find('hemi-R') != -1: hemi = 'hemi-R'
    else: hemi = None


    # prepoc files name

    for task in tasks : 

        glm_files_task = [file for file in preproc_fn_runs if task in file]
        # Contrast
        if task == 'SacLoc':
            cond1_label, cond2_label = ['Sac'], ['Fix']
        elif task == 'PurLoc':
            cond1_label, cond2_label = ['Pur'], ['Fix']
            

    
    
    
        events_tot = pd.DataFrame()
        confounds_tot = pd.DataFrame()
        preproc_data_tot = []
        n_runs = 0
        for preproc_fn_idx, preproc_fn in enumerate(glm_files_task) :
            match = re.search(r'_run-(\d+)_', preproc_fn)
            run_num = 'run-{}'.format(match.group(1))
            n_runs += 1

        
        
        

        
            
            # find the events and confounds files 
            event_dir = '{}/{}/{}/{}/func/'.format(main_dir, 
                                                   project_dir, 
                                                   subject, 
                                                   func_session)
            
            con_dir = '{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func'.format(main_dir, 
                                                                       project_dir, 
                                                                       subject, 
                                                                       func_session)
            
                
        

            # Find the event files 
            event_file = glob.glob("{}/{}_{}_task-{}_{}_events.tsv".format(event_dir, 
                                                                              subject, 
                                                                              func_session, 
                                                                              task,
                                                                              run_num))
            # Finf the confounds files 
            con_file = glob.glob('{}/{}_{}_task-{}_{}_desc-confounds_timeseries.tsv'.format(con_dir, 
                                                                                 subject, 
                                                                                 func_session, 
                                                                                 task,
                                                                                 run_num))
            
            deb()
            
        
            confounds = pd.read_table(con_file[0])[confounds_list].dropna(axis=1)
            confounds_tot = pd.concat([confounds_tot, confounds], ignore_index=True)
        
            # make the designe matrixe  
            events = eventsMatrix(design_file=event_file[0], task=task, tr=TR)
            events_tot = pd.concat([events_tot, events], ignore_index=True)
      
            # Load data
            print('load',preproc_fn)
            preproc_img, preproc_data = load_surface(fn=preproc_fn)
            preproc_data_tot.append(preproc_data)
            

        preproc_data_tot = np.concatenate(preproc_data_tot, axis=0)

        
        
        
        events_tot = events_tot.sort_values(by='onset').reset_index(drop=True)
        events_tot['onset_shifted'] = 0
        for i in range(1, len(events_tot)):
            events_tot.loc[i, 'onset_shifted'] = events_tot.loc[i - 1, 'onset_shifted'] + events_tot.loc[i - 1, 'duration']
        
        
        events_tot = events_tot.drop(columns='onset').rename(columns={'onset_shifted': 'onset'})
      

            
            
        
        frame_times = np.arange(TRs*n_runs) * TR
        design_matrix = make_first_level_design_matrix(frame_times,
                                                   events=events_tot,
                                                   hrf_model='spm',
                                                   drift_model=None,
                                                   add_regs=confounds_tot)

    







        # Run the glm
        labels, estimates = run_glm(preproc_data_tot, design_matrix.values, noise_model="ols")
        
        # extract glm predictions and r2       
        glm_pred, glm_r2 = extract_predictions_r2 (labels=labels,
                                                   estimate=estimates,
                                                   source_data=preproc_data_tot)
    
        
        # Compute the contrasts 
        for contrast_num, contrast in enumerate(zip(cond1_label,cond2_label)):
            # make contrasts
            contrast_values = (design_matrix.columns == contrast[0]) * 1.0 -(design_matrix.columns == contrast[1])
            # compute contrasts
            eff = compute_contrast(labels, estimates, contrast_values,contrast_type='t')
        
            
            # compute the derivatives               
            z_map = eff.z_score()
            z_p_map = 2*(1 - stats.norm.cdf(abs(z_map)))
            fdr_th = fdr_threshold(z_map, glm_alpha)
            fdr = z_map
            fdr *= (z_map > fdr_th)
            fdr_p_map = 2*(1 - stats.norm.cdf(abs(fdr)))
            
            if contrast_num:
                fit = np.vstack((fit,z_map,z_p_map,fdr,fdr_p_map))
            else:                 
                fit = np.vstack((z_map,z_p_map,fdr,fdr_p_map,glm_r2))
                

           
            

        if hemi:
            glm_fit_fn = "{}/{}/derivatives/pp_data/{}/fsnative/glm/glm_fit/{}_task-{}_{}_space-fsnative_dct_glm-fit_concat.func.gii".format(
                main_dir, project_dir, subject, subject, task ,hemi)
            glm_pred_fn = "{}/{}/derivatives/pp_data/{}/fsnative/glm/glm_fit/{}_task-{}_{}_space-fsnative_dct_glm-pred_concat.func.gii".format(
                main_dir, project_dir, subject, subject, task ,hemi)
   
            
        else:
            glm_fit_fn = "{}/{}/derivatives/pp_data/{}/170k/glm/glm_fit/{}_task-{}_space-fsLR_den-170k_dct_glm-fit_concat.dtseries.nii".format(
                main_dir, project_dir, subject, subject,task)
            glm_pred_fn = "{}/{}/derivatives/pp_data/{}/170k/glm/glm_fit/{}_task-{}_space-fsLR_den-170k_dct_glm-pred_concat.dtseries.nii".format(
                main_dir, project_dir, subject, subject,task)
        # export fit
        maps_names = ['z_map','z_p_map','fdr','fdr_p_map','rsquare_map']

        fit_img = make_surface_image(data=fit, 
                                     source_img=preproc_img, 
                                     maps_names=maps_names)

        
        nb.save(fit_img,glm_fit_fn) 
        
        # export pred
        pred_img = make_surface_image(data=glm_pred, 
                                      source_img=preproc_img)

        nb.save(pred_img,glm_pred_fn) 
        
        print('{} is done'.format(preproc_fn))



        
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))

        


