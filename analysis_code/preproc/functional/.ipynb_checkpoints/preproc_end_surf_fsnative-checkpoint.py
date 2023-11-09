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
>> cd /home/mszinte/projects/stereo_prf/analysis_code/preproc/functional/
2. run python command
python preproc_end.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python preproc_end_surf.py /scratch/mszinte/data retinoMaps sub-01 327
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
import ipdb
import numpy as np
import nibabel as nb
import itertools as it
from scipy import stats
from nilearn import signal

from nilearn.glm.first_level.design_matrix import _cosine_drift
trans_cmd = 'rsync -avuz --progress'
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
high_pass_threshold = analysis_info['high_pass_threshold'] 
high_pass_type = analysis_info['high_pass_type'] 
sessions = analysis_info['session']

# main_dir = '/Users/uriel/disks/meso_shared'
# project_dir = 'RetinoMaps'
# subject = 'sub-02'
# group = '327'
# TR = 1.3
# tasks = 'pRF'
# high_pass_threshold = 0.01
# high_pass_type = 'dct'
# sessions = ['ses-01']

# # load settings
# with open('/Users/uriel/disks/meso_H/projects/RetinoMaps/analysis_code/settings.json') as f:
#     json_s = f.read()
#     analysis_info = json.loads(json_s)
# TR = analysis_info['TR']

# high_pass_threshold = analysis_info['high_pass_threshold'] 
# high_pass_type = analysis_info['high_pass_type'] 



for session in sessions :
    

    
    if session == 'ses-01':
        tasks = ['pRF']
    else : 
        tasks = ["rest","pMF","PurLoc","SacLoc","PurVELoc","SacVELoc"]
         
    fmriprep_dir = "{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/".format(main_dir, project_dir, subject, session)    
    fmriprep_func_RH_fns = glob.glob("{}/*_hemi-R_space-fsnative_bold.func.gii".format(fmriprep_dir))
    fmriprep_func_LH_fns = glob.glob("{}/*_hemi-L_space-fsnative_bold.func.gii".format(fmriprep_dir))
    pp_data_func_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/fsnative".format(main_dir, project_dir, subject)
    os.makedirs(pp_data_func_dir, exist_ok=True)
    
    # High pass filtering and z-scoring
    print("high-pass filtering...")
    for func_fn_R,func_fn_L in zip(fmriprep_func_RH_fns,fmriprep_func_LH_fns) :
        # Load data
        surface_data_R_img = nb.load(func_fn_R)
        surface_data_R_header = surface_data_R_img.header
        surface_data_R_meta = surface_data_R_img.meta
        
        surface_data_L_img = nb.load(func_fn_L)
        surface_data_L_header = surface_data_L_img.header
        surface_data_L_meta = surface_data_L_img.meta
        
        
        # Right hemisphere
        surf_data_R = [x.data for x in surface_data_R_img.darrays]
        surf_data_R = np.vstack(surf_data_R)
        
        # Left hemisphere    
        surf_data_L = [x.data for x in surface_data_L_img.darrays]
        surf_data_L = np.vstack(surf_data_L)
        
        # High pass filtering Left hemisphere
        n_vol_L = surf_data_L.shape[0]
        ft_L = np.linspace(0.5 * TR, (n_vol_L + 0.5) * TR, n_vol_L, endpoint=False)
        hp_set_L = _cosine_drift(high_pass_threshold, ft_L)
        surf_data_L = signal.clean(surf_data_L, detrend=False, standardize=True, confounds=hp_set_L)
        
        # High pass filtering Right hemisphere
        n_vol_R = surf_data_R.shape[0]
        ft_R = np.linspace(0.5 * TR, (n_vol_R + 0.5) * TR, n_vol_R, endpoint=False)
        hp_set_R = _cosine_drift(high_pass_threshold, ft_R)
        surf_data_R = signal.clean(surf_data_R, detrend=False, standardize=True, confounds=hp_set_R)
    
        #Export left hemisphere filtering data
        flt_im_L = nb.gifti.GiftiImage(header = surface_data_L_header, meta = surface_data_L_meta)
        for flt_data_L in surf_data_L:
            flt_darray_L = nb.gifti.GiftiDataArray(flt_data_L)
            flt_im_L.add_gifti_data_array(flt_darray_L)
            
        out_flt_file_L = "{}/{}_{}.func.gii".format(pp_data_func_dir,func_fn_L.split('/')[-1][:-9],high_pass_type) 
        nb.save(flt_im_L, out_flt_file_L)
        
        
        #Export right hemisphere filtering data
        flt_im_R = nb.gifti.GiftiImage(header = surface_data_R_header, meta = surface_data_R_meta)
        for flt_data_R in surf_data_R:
            flt_darray_R = nb.gifti.GiftiDataArray(flt_data_R)
            flt_im_R.add_gifti_data_array(flt_darray_R)
        
        out_flt_file_R = "{}/{}_{}.func.gii".format(pp_data_func_dir,func_fn_R.split('/')[-1][:-9],high_pass_type)   
        nb.save(flt_im_R, out_flt_file_R)
 
        
 
    for task in tasks:
        print(task)

        # preproc tasks runs
        preproc_files_R = glob.glob("{}/*_task-{}_*_hemi-R_space-fsnative_bold_{}.func.gii".format(pp_data_func_dir,task, high_pass_type))
        preproc_files_L = glob.glob("{}/*_task-{}_*_hemi-L_space-fsnative_bold_{}.func.gii".format(pp_data_func_dir,task, high_pass_type))
        
        
        # make correlation directory
        corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/fsnative".format(main_dir, project_dir, subject)
        os.makedirs(corr_dir, exist_ok=True)

        # output files 
        cor_file_R = "{}/{}_task-{}_hemi-R_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,high_pass_type)
        cor_file_L = "{}/{}_task-{}_hemi-L_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,high_pass_type)
        
        
        # load preproc files to have meta and header
        preproc_img_L = nb.load(preproc_files_L[0])
        preproc_data_L = [x.data for x in preproc_img_L.darrays]
        preproc_data_L = np.vstack(preproc_data_L)
        
        preproc_img_header_L = preproc_img_L.header
        preproc_img_meta_L = preproc_img_L.meta
        
        
        
        # make empty right image
        preproc_img_R = nb.load(preproc_files_R[0])
        preproc_data_R = [x.data for x in preproc_img_R.darrays]
        preproc_data_R = np.vstack(preproc_data_R)
        
        preproc_img_header_R = preproc_img_R.header
        preproc_img_meta_R = preproc_img_R.meta
       


        

        # Correlations
        print('starting correlations')
        
        # compute the combination 
        combis_L = list(it.combinations(preproc_files_L, 2))
        combis_R = list(it.combinations(preproc_files_R, 2))
        
        # Left hemisphere
        cor_final_L = np.zeros((preproc_data_L.shape[1]))
        for combi_L in combis_L:
            task_cor_L = np.zeros((preproc_data_L.shape[1]))
            
            a_img_L = nb.load(combi_L[0])
            a_data_L = [x.data for x in a_img_L.darrays]
            a_data_L = np.vstack(a_data_L)
    
            b_img_L = nb.load(combi_L[1])
            b_data_L = [x.data for x in b_img_L.darrays]
            b_data_L = np.vstack(b_data_L)
            
            for vertice_L in range(a_data_L.shape[1]):
                corr_L, _ = stats.pearsonr(a_data_L[:, vertice_L], b_data_L[:, vertice_L])
                task_cor_L[vertice_L] = corr_L
        
            cor_final_L += task_cor_L / len(combi_L)
        
        print('starting exportation')
        
        # export left hemisphere correlations file
        corr_img_L = nb.gifti.GiftiImage(header=preproc_img_header_L, meta=preproc_img_meta_L)
        for corr_data_L in cor_final_L:
            corr_darray_L = nb.gifti.GiftiDataArray(corr_data_L)
            corr_img_L.add_gifti_data_array(corr_darray_L)
        cor_file_L = "{}/{}_task-{}_hemi-L_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task, high_pass_type)
        nb.save(corr_img_L, cor_file_L)
        
        print('left hemisphere correlation done')
        
        # Right hemisphere
        cor_final_R = np.zeros((preproc_data_R.shape[1]))
        for combi_R in combis_R:
            task_cor_R = np.zeros((preproc_data_R.shape[1]))
            
            a_img_R = nb.load(combi_R[0])
            a_data_R = [x.data for x in a_img_R.darrays]
            a_data_R = np.vstack(a_data_R)
    
            b_img_R = nb.load(combi_R[1])
            b_data_R = [x.data for x in b_img_R.darrays]
            b_data_R = np.vstack(b_data_R)
        
            for vertice_R in range(a_data_R.shape[1]):

        
                corr_R, _ = stats.pearsonr(a_data_R[:, vertice_R], b_data_R[:, vertice_R])
                task_cor_R[vertice_R] = corr_R
        
            cor_final_R += task_cor_R / len(combi_R)
        
        # export right hemisphere correlations file
        corr_img_R = nb.gifti.GiftiImage(header=preproc_img_header_R, meta=preproc_img_meta_R)
        for corr_data_R in cor_final_R:
            corr_darray_R = nb.gifti.GiftiDataArray(corr_data_R)
            corr_img_R.add_gifti_data_array(corr_darray_R)
        cor_file_R = "{}/{}_task-{}_hemi-R_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task, high_pass_type)
        nb.save(corr_img_R, cor_file_R)
        print('right hemisphere correlation done')


        # Averaging 
        print('averaging')
        # make avg dir
        avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/fsnative".format(main_dir, project_dir, subject)
        os.makedirs(avg_dir, exist_ok=True)
        
        
        # avg files 
        avg_file_R = "{}/{}_task-{}_hemi-R_fmriprep_{}_avg_bold.func.gii".format(avg_dir, subject, task,high_pass_type)
        avg_file_L = "{}/{}_task-{}_hemi-L_fmriprep_{}_avg_bold.func.gii".format(avg_dir, subject, task,high_pass_type)
        
        

         
        print("averaging...")
        for preproc_file_L, preproc_file_R in zip(preproc_files_L,preproc_files_R):
            
            data_avg_R = np.zeros(preproc_data_R.shape)
            data_avg_L = np.zeros(preproc_data_L.shape) 

            print('add: {}'.format(preproc_file_L))
            data_val_L = []
            data_val_R = []
            
            # Load left hemisphere data
            data_val_img_L = nb.load(preproc_file_L)
            
            
            header_val_img_L = data_val_img_L.header
            meta_val_img_L = data_val_img_L.meta
            
            data_val_L = [x.data for x in data_val_img_L.darrays]
            data_val_L = np.vstack(data_val_L)
            
            # Load right hemisphere data
            data_val_img_R = nb.load(preproc_file_R)
            
            header_val_img_R = data_val_img_R.header
            meta_val_img_R = data_val_img_R.meta
            
            data_val_R = [x.data for x in data_val_img_R.darrays]
            data_val_R = np.vstack(data_val_R)
            
            # Averaging 
            data_avg_L += data_val_L/len(preproc_files_L)
            data_avg_R += data_val_R/len(preproc_files_R)
            
            # export left hesmisphere averaging data
            avg_img_L = nb.gifti.GiftiImage(header = header_val_img_L, meta= meta_val_img_L)
            for avg_data_L in data_avg_L:
                avg_darray_L = nb.gifti.GiftiDataArray(avg_data_L)
                avg_img_L.add_gifti_data_array(avg_darray_L)  
            nb.save(avg_img_L, avg_file_L)
            
            # export right hesmisphere averaging data
            avg_img_R = nb.gifti.GiftiImage(header = header_val_img_R, meta= meta_val_img_R)
            for avg_data_R in data_avg_R:
                avg_darray_R = nb.gifti.GiftiDataArray(avg_data_R)
                avg_img_L.add_gifti_data_array(avg_darray_R)  
            nb.save(avg_img_R, avg_file_R)
                
            # except :      
            #     print("An error occure during averaging task {task}".format(task = task))
        
            
        # Leave-one-out averages
        if len(preproc_files_L):
            combi_L = list(it.combinations(preproc_files_L, len(preproc_files_L)-1))
        
        if len(preproc_files_R):
            combi_R = list(it.combinations(preproc_files_R, len(preproc_files_R)-1))
        
        
        # Left hemisphere 
        for loo_num, avg_runs in enumerate(combi_L):
            loo_avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_loo_avg/fsnative".format(main_dir, project_dir, subject)
            os.makedirs(loo_avg_dir, exist_ok=True)
            
            
            # try:
            print(loo_num)


            print("loo_avg-{}".format(loo_num+1))
            
            # compute average between loo runs
            loo_avg_file_L = "{}/{}_task-{}_hemi-L_fmriprep_{}_avg_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)


        
        
            preproc_val_L = nb.load(preproc_files_L[0])
            preproc_data_L = [x.data for x in preproc_val_L.darrays]
            preproc_data_L = np.vstack(preproc_data_L)
            data_loo_avg_L = np.zeros(preproc_data_L.shape)
        
            for avg_run in avg_runs:
                print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
                data_val_L = []
                data_val_img_L = nb.load(avg_run)
                data_val_L = [x.data for x in data_val_img_L.darrays]
                data_val_L = np.vstack(data_val_L)
                data_loo_avg_L += data_val_L/len(avg_runs)
                
            
            loo_avg_img_L = nb.gifti.GiftiImage(header = header_val_img_L, meta= meta_val_img_L)
            for data_loo_L in data_loo_avg_L:
                darray_loo_L = nb.gifti.GiftiDataArray(data_loo_L)
                loo_avg_img_L.add_gifti_data_array(darray_loo_L)
            
                                      
            nb.save(loo_avg_img_L, loo_avg_file_L)        
        
            # copy loo run (left one out run)
            for loo in preproc_files_L:
                if loo not in avg_runs:
                    loo_file_L =  "{}/{}_task-{}_hemi-L_fmriprep_{}_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)
                    print("loo: {}".format(loo))
                    os.system("{} {} {}".format(trans_cmd, loo, loo_file_L))

            # except : 
            #     print("An error occure during loo task {task} left hemisphere".format(task = task))
        
        # Right hemisphere 
        for loo_num, avg_runs in enumerate(combi_R):
            # try : 
            print(loo_num)


            print("loo_avg-{}".format(loo_num+1))
            
            # compute average between loo runs
            loo_avg_file_R = "{}/{}_task-{}_hemi-R_fmriprep_{}_avg_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)


        
        
            preproc_val_R = nb.load(preproc_files_R[0])
            preproc_data_R = [x.data for x in preproc_val_R.darrays]
            preproc_data_R = np.vstack(preproc_data_R)
            data_loo_avg_R = np.zeros(preproc_data_R.shape)
        
            for avg_run in avg_runs:
                print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
                data_val_R = []
                data_val_img_R = nb.load(avg_run)
                data_val_R = [x.data for x in data_val_img_R.darrays]
                data_val_R = np.vstack(data_val_R)
                data_loo_avg_R += data_val_R/len(avg_runs)
                
            
            loo_avg_img_R = nb.gifti.GiftiImage(header = header_val_img_R, meta= meta_val_img_R)
            for data_loo_R in data_loo_avg_R:
                darray_loo_R = nb.gifti.GiftiDataArray(data_loo_R)
                loo_avg_img_R.add_gifti_data_array(darray_loo_R)
            
                                      
            nb.save(loo_avg_img_R, loo_avg_file_R)        
        
            # copy loo run (left one out run)
            for loo in preproc_files_R:
                if loo not in avg_runs:
                    loo_file_R =  "{}/{}_task-{}_hemi-R_fmriprep_{}_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)
                    print("loo: {}".format(loo))
                    os.system("{} {} {}".format(trans_cmd, loo, loo_file_L))


