"""
-----------------------------------------------------------------------------------------
preproc_end_surf.py
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
>> cd ~/projects/RetinoMaps/analysis_code/preproc/functional/
2. run python command
python preproc_end_surf.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python preproc_end_surf.py /scratch/mszinte/data RetinoMaps sub-02 327
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
import shutil
from nilearn.glm.first_level.design_matrix import _cosine_drift
import datetime



sys.path.append("{}/../../utils".format(os.getcwd()))

sys.path.append('/Users/uriel/disks/meso_H/projects/RetinoMaps/analysis_code/utils')
from surface_utils import load_surface , make_surface_image
deb = ipdb.set_trace
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
high_pass_threshold = analysis_info['high_pass_threshold'] 
high_pass_type = analysis_info['high_pass_type'] 
sessions = analysis_info['session']
tasks = analysis_info['task_names']



formats = ['fsnative','170k']
extensions = ['func.gii','dtseries.nii']


for format_, extension in zip(formats,extensions):
    # make  directorys
    flt_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct".format(
        main_dir, project_dir, subject)
    os.makedirs('{}/{}'.format(flt_dir,format_), exist_ok=True)
    
    corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/{}".format(
        main_dir, project_dir, subject,format_)
    os.makedirs(corr_dir, exist_ok=True)
    
    avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/{}".format(
        main_dir, project_dir, subject,format_)
    os.makedirs(avg_dir, exist_ok=True)
    
    loo_avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_loo_avg/{}".format(
        main_dir, project_dir, subject,format_)
    os.makedirs(loo_avg_dir, exist_ok=True)


   # # DCT correction
   # # -----------------------
   #  for session in sessions :

   #      fmriprep_dir = "{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/".format(
   #          main_dir, project_dir, subject, session)
   #      fmriprep_func_fns = glob.glob("{}/*-{}*.{}".format(
   #          fmriprep_dir,format_,extension)) 

    
                
    

   #      for func_fn in fmriprep_func_fns :
   #          # make output filtered files
   #          flt_data_fn = func_fn.split('/')[-1]
   #          flt_data_fn =flt_data_fn.replace('bold','dct_bold')


   #          # Load data
   #          surf_img, surf_data = load_surface(fn=func_fn)
            
   #          # High pass filtering 
   #          n_tr = surf_data.shape[0]
   #          ft = np.linspace(0.5 * TR, (n_tr + 0.5) * TR, n_tr, endpoint=False)
   #          hp_set = _cosine_drift(high_pass_threshold, ft)
   #          surf_data = signal.clean(surf_data, detrend=False,
   #                                   standardize=False, confounds=hp_set)
            
   #          # Compute the Z-score 
   #          surf_data =  (surf_data - np.mean(surf_data, axis=0)) / np.std(surf_data, axis=0)
        
   #          # Make an image with the preproceced data
            
   #          flt_img = make_surface_image(data=surf_data,source_img=surf_img)   
   #          nb.save(flt_img, '{}/{}/{}'.format(flt_dir,format_,flt_data_fn))
     
    
    
# find all the filtered files 
preproc_files_tot = glob.glob("{}/{}/*_*.{}".format(
    flt_dir,formats[0], extensions[0])) + glob.glob("{}/{}/*_*.{}".format(
        flt_dir,formats[1], extensions[1]))
        
        
# split filtered files  deppending of their nature 
preproc_fsnative_hemi_L = []
preproc_fsnative_hemi_R = []
preproc_170k = []

for subtype in preproc_files_tot:
    if "hemi-L" in subtype:
        preproc_fsnative_hemi_L.append(subtype)
        
    elif "hemi-R" in subtype:
        preproc_fsnative_hemi_R.append(subtype)
        
    elif "170k" in subtype:
        preproc_170k.append(subtype)
        
preproc_files_list = [preproc_fsnative_hemi_L,preproc_fsnative_hemi_R,preproc_170k]


# run correlations , averagin and loo averagin for each subtype 
for preproc_files in preproc_files_list:
    for task in tasks:

        

        # defind output files names 
        if preproc_files[0].find('hemi-L') != -1:
            hemi = 'hemi-L'
        elif preproc_files[0].find('hemi-R') != -1:
            hemi = 'hemi-R'
        else:
            hemi = None  
            

        if hemi:
            cor_file = "{}/{}_task-{}_{}_fmriprep_{}_correlations_bold.func.gii".format(
                corr_dir, subject, task, hemi, high_pass_type)

            avg_file = "{}/{}_task-{}_{}_fmriprep_{}_avg_bold.func.gii".format(
                avg_dir, subject, task, hemi, high_pass_type)

            loo_avg_file = "{}/{}_task-{}_{}_fmriprep_{}_avg_loo-{}_bold.func.gii".format(
                loo_avg_dir, subject, task, hemi, high_pass_type, 'loo_num')

            loo_file = "{}/{}_task-{}_{}_fmriprep_{}_loo-{}_bold.func.gii".format(
                loo_avg_dir, subject, task, hemi, high_pass_type, 'loo_num')
        else:
            cor_file = "{}/{}_task-{}_fmriprep_{}_correlations_bold.dtseries.nii".format(
                corr_dir, subject, task, high_pass_type)

            avg_file = "{}/{}_task-{}_fmriprep_{}_avg_bold.dtseries.nii".format(
                avg_dir, subject, task, high_pass_type)

            loo_avg_file = "{}/{}_task-{}_fmriprep_{}_avg_loo-{}_bold.dtseries.nii".format(
                loo_avg_dir, subject, task, high_pass_type, 'loo_num')

            loo_file = "{}/{}_task-{}_fmriprep_{}_loo-{}_bold.dtseries.nii".format(
                loo_avg_dir, subject, task, high_pass_type, 'loo_num')

        # load preproc files to have meta and header
        preproc_img, preproc_data = load_surface(fn=preproc_files[0])
        
    
        # Correlation computation
        # -----------------------
        print('starting correlations')
        
        # compute the combination 
        combis = list(it.combinations(preproc_files, 2))
        
        #  load data and comute the correaltions
        cor_final = np.zeros((1,preproc_data.shape[1]))
        for combi in combis:
            task_cor = np.zeros((preproc_data.shape[1]))
    
            a_img, a_data = load_surface(fn=combi[0])
            b_img, b_data = load_surface(fn=combi[1])
            
            for vertice in range(a_data.shape[1]):
                corr, _ = stats.pearsonr(a_data[:, vertice], 
                                              b_data[:, vertice])
                task_cor[vertice] = corr
        
            cor_final += task_cor / len(combis)
        

    
        corr_img = make_surface_image(data=cor_final,source_img=preproc_img)   
        nb.save(corr_img, cor_file)
        
        
        # Averaging computation
        # ---------------------
        print('averaging')
        
    
        data_avg = np.zeros(preproc_data.shape)
        
        for preproc_file in preproc_files:
            
            # Load data
            preproc_img, preproc_data = load_surface(fn=preproc_file)
            
            # Averaging 
            data_avg += preproc_data/len(preproc_files)
    
        # export averaging data
        avg_img = make_surface_image(data =data_avg, source_img = preproc_img)
        nb.save(avg_img, avg_file)
            
        
        # Leave-one-out averages computation
        # ---------------------
        
        if len(preproc_files):
            combi = list(it.combinations(preproc_files, len(preproc_files)-1))
        
        for loo_num, avg_runs in enumerate(combi):
    
            
            print("loo_avg-{}".format(loo_num+1))
            
            # compute average between loo runs
            loo_avg_file = loo_avg_file.replace('loo_num','{}'.format(loo_num+1))
    
            # load data and make the loo_avg object
            preproc_img, preproc_data = load_surface(fn=preproc_files[0])
            data_loo_avg = np.zeros(preproc_data.shape)
        
            # compute leave on out averagin
            for avg_run in avg_runs:
                print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
                
                preproc_img, preproc_data = load_surface(fn=avg_run)
                
    
                data_loo_avg += preproc_data/len(avg_runs)
                
            # export leave one out file 
            loo_avg_img = make_surface_image(data = data_loo_avg, source_img=preproc_img)
            nb.save(loo_avg_img, loo_avg_file)        
        
            # copy loo run (left one out run)
            for loo in preproc_files:
                if loo not in avg_runs:
                    loo_file = loo_file.replace('loo_num','{}'.format(loo_num+1))
    
                    print("loo: {}".format(loo))
                    shutil.copyfile(loo, loo_file)
                        
# rename correlations maps for fsnative
fsnative_corr_dir  = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/fsnative".format(
    main_dir, project_dir, subject)
cors_fn = glob.glob('{}/*.func.gii'.format(fsnative_corr_dir))

for cor_fn in cors_fn :
    os.system('wb_command -set-map-names {} -map {} {}'.format(
        cor_fn,1,'runs_correlations'))

# Anatomy
print("getting anatomy...")
orig_dir_anat = "{}/{}/derivatives/fmriprep/fmriprep/{}/ses-01/anat/".format(
    main_dir, project_dir, subject)
pycortex_flat_dir = '{}/{}/derivatives/pp_data/cortex/db/{}/surfaces'.format(
    main_dir,project_dir,subject)
anat_files = glob.glob("{}/*.surf.gii".format(orig_dir_anat))

dest_dir_anat = "{}/{}/derivatives/pp_data/{}/anat".format(
    main_dir, project_dir, subject)
os.makedirs(dest_dir_anat, exist_ok=True)
hemis = ['L','R']
# load flat data and change medadata to make them readable by wb_view
for hemi in hemis : 
    if hemi == 'L' :
        flat_img_l = nb.load('{}/flat_lh.gii'.format(pycortex_flat_dir))
        flat_img_l.darrays[0].meta['AnatomicalStructurePrimary'] = 'CortexLeft'
        flat_img_l.darrays[0].meta['GeometricType']= 'Flat'
        nb.save(flat_img_l,'{}/{}_flat_lh.surf.gii'.format(
            dest_dir_anat,subject))
        
    elif hemi == 'R' :
        flat_img_r = nb.load('{}/flat_rh.gii'.format(pycortex_flat_dir))
        flat_img_r.darrays[0].meta['AnatomicalStructurePrimary'] = 'CortexRight'
        flat_img_r.darrays[0].meta['GeometricType']= 'Flat'
        nb.save(flat_img_r,'{}/{}_flat_rh.surf.gii'.format(
            dest_dir_anat,subject))

# import surface anat data 
for orig_file in anat_files:
    file_name = os.path.basename(orig_file)
    dest_file = os.path.join(dest_dir_anat, file_name)
    shutil.copyfile(orig_file, dest_file)
    
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))
    
    
        
                                                                   
