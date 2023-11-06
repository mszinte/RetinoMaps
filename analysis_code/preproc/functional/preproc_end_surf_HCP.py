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
import ipdb
import platform
import numpy as np
import nibabel as nb
import itertools as it
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
# subject = 'sub-02b'
# group = '327'
# session = 'ses-01'


# file_path = '/Users/uriel/disks/meso_shared/RetinoMaps/derivatives/fmriprep/fmriprep/sub-02b/ses-02/func'
# main_dir = '/Users/uriel/disks/meso_shared/'
# project_dir = 'RetinoMaps'
# file_name ='sub-02_ses-02_task-pMF_run-01_hemi-L_space-fsaverage_bold.func.gii'
# subject = 'sub-03'

# # load settings
# with open('/Users/uriel/disks/meso_H/projects/RetinoMaps/analysis_code/settings.json') as f:
#     json_s = f.read()
#     analysis_info = json.loads(json_s)
# TR = analysis_info['TR']
# #tasks = analysis_info['task_names']
# high_pass_threshold = analysis_info['high_pass_threshold'] 
# high_pass_type = analysis_info['high_pass_type'] 
# sessions = analysis_info['session']



for session in sessions : 
    
    if session == 'ses-01':
        tasks = ['pRF']
    else : 
        tasks = ["rest","pMF","PurLoc","SacLoc","PurVELoc","SacVELoc"]
    
    
    # Get fmriprep filenames
    fmriprep_dir = "{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/".format(main_dir, project_dir, subject, session)
    fmriprep_func_fns = glob.glob("{}/*_space-fsLR_den-170k_bold.dtseries.nii".format(fmriprep_dir))
    
    pp_data_func_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/HCP_170k".format(main_dir, project_dir, subject)
    os.makedirs(pp_data_func_dir, exist_ok=True)
    
    # High pass filtering and z-scoring
    print("high-pass filtering...")
    for func_fn in fmriprep_func_fns:
        surf_im = nb.load(func_fn)
        surf_data = surf_im.get_fdata()
    
    
        
    
        n_vol = surf_data.shape[0]
        ft = np.linspace(0.5 * TR, (n_vol + 0.5) * TR, n_vol, endpoint=False)
        hp_set = _cosine_drift(high_pass_threshold, ft)
        flt_surf_data = signal.clean(surf_data, detrend=False, standardize=True, confounds=hp_set)
        
    
        flt_surf_im = nb.cifti2.cifti2.Cifti2Image(dataobj=flt_surf_data, header=surf_im.header)
        out_flt_file = "{}/{}_{}.dtseries.nii".format(pp_data_func_dir,func_fn.split('/')[-1][:-13],high_pass_type) 
        nb.save(flt_surf_im, out_flt_file)
        
        
        
    for task in tasks:    
        # Average tasks runs
        preproc_files = glob.glob("{}/*_task-{}_*_space-fsLR_den-170k_bold_{}.dtseries.nii".format(pp_data_func_dir,task, high_pass_type))
        print(preproc_files)
    
    
        avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/HCP_170k".format(main_dir, project_dir, subject)
        os.makedirs(avg_dir, exist_ok=True)
        
        avg_file = "{}/{}_task-{}_space-fsLR_den-170k_bold_{}_avg.dtseries.nii".format(avg_dir, subject, task,high_pass_type)
        img = nb.load(preproc_files[0])
        data_avg = np.zeros(img.shape)
        
        print("averaging...")
        for file in preproc_files:
            print('add: {}'.format(file))
            data_val = []
            data_val_img = nb.load(file)
            data_val = data_val_img.get_fdata()
            data_avg += data_val/len(preproc_files)
        
    
        
        avg_img = nb.cifti2.cifti2.Cifti2Image(dataobj=data_avg, header=img.header)
        nb.save(avg_img, avg_file)
        
        
        # Leave-one-out averages
        if len(preproc_files):
            combi = list(it.combinations(preproc_files, len(preproc_files)-1))
       
        loo_avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_loo_avg/HCP_170k".format(main_dir, project_dir, subject)
        os.makedirs(loo_avg_dir, exist_ok=True)

        
        for loo_num, avg_runs in enumerate(combi):
            print("loo_avg-{}".format(loo_num+1))
        
            # compute average between loo runs
            loo_avg_file = "{}/{}_task-{}_fmriprep_bold_{}_avg_loo-{}.dtseries.nii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)
            
            img = nb.load(preproc_files[0])
            data_loo_avg = np.zeros(img.shape)
        
            for avg_run in avg_runs:
                print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
                data_val = []
                data_val_img = nb.load(avg_run)
                data_val = data_val_img.get_fdata()
                data_loo_avg += data_val/len(avg_runs)
        
            
            loo_avg_img = nb.cifti2.cifti2.Cifti2Image(dataobj=data_loo_avg, header=img.header)
            nb.save(loo_avg_img, loo_avg_file)
        
            # copy loo run (left one out run)
            for loo in preproc_files:
                if loo not in avg_runs:
                    loo_file = "{}/{}_task-{}_fmriprep_bold_{}_loo-{}.dtseries.nii".format(loo_avg_dir, subject,task,high_pass_type, loo_num+1)
                    print("loo: {}".format(loo))
                    os.system("{} {} {}".format(trans_cmd, loo, loo_file))
                                                    
    # Anatomy
    print("getting anatomy...")
    output_files = ['dseg','desc-preproc_T1w','desc-aparcaseg_dseg','desc-aseg_dseg','desc-brain_mask']
    orig_dir_anat = "{}/{}/derivatives/fmriprep/fmriprep/{}/ses-01/anat/".format(main_dir, project_dir, subject)
    dest_dir_anat = "{}/{}/derivatives/pp_data/{}/anat".format(main_dir, project_dir, subject)
    os.makedirs(dest_dir_anat,exist_ok=True)
    
    for output_file in output_files:
        orig_file = "{}/{}_{}_{}.nii.gz".format(orig_dir_anat, subject, session, output_file)
        dest_file = "{}/{}_{}.nii.gz".format(dest_dir_anat, subject, output_file)
        os.system("{} {} {}".format(trans_cmd, orig_file, dest_file))
        

