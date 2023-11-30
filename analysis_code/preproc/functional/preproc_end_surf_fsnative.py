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
>> cd ~/projects/RetinoMaps/analysis_code/preproc/functional/
2. run python command
python preproc_end.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python preproc_end_surf_fsnative.py /scratch/mszinte/data RetinoMaps sub-02 327
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
trans_cmd = 'rsync -avuz --progress'

sys.path.append("{}/../../utils".format(os.getcwd()))
from gifti_utils import make_giti_image
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
hemis = ['L','R']



for session in sessions :
    
    if session == 'ses-01':
        tasks = ['pRF']
    else : 
        tasks = ["rest","pMF","PurLoc","SacLoc","PurVELoc","SacVELoc"]
    
   
    
    fmriprep_dir = "{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/".format(main_dir, project_dir, subject, session)
    pp_data_func_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/fsnative".format(main_dir, project_dir, subject)
    os.makedirs(pp_data_func_dir, exist_ok=True)
    
    for hemi in hemis :
        fmriprep_func_hemi_fns = glob.glob("{}/*_hemi-{}_space-fsnative_bold.func.gii".format(fmriprep_dir,hemi))

        
        # High pass filtering and z-scoring
        print("high-pass filtering...")
        for func_fn_hemi in fmriprep_func_hemi_fns :
            

            # Load data
            surface_data_hemi_img = nb.load(func_fn_hemi)

            
 
            surf_data_hemi = [x.data for x in surface_data_hemi_img.darrays]
            surf_data_hemi = np.vstack(surf_data_hemi)
            
            # High pass filtering Left hemisphere
            n_vol_hemi = surf_data_hemi.shape[0]
            ft_hemi = np.linspace(0.5 * TR, (n_vol_hemi + 0.5) * TR, n_vol_hemi, endpoint=False)
            hp_set_hemi = _cosine_drift(high_pass_threshold, ft_hemi)
            surf_data_hemi = signal.clean(surf_data_hemi, detrend=False, standardize=False, confounds=hp_set_hemi)
            
            # Compute the Z-score 
            surf_data_hemi =  (surf_data_hemi - np.mean(surf_data_hemi, axis=0)) / np.std(surf_data_hemi, axis=0)
            
        
            # Make a GIFTI image with the preproceced data
            flt_im_hemi = make_giti_image(surface_data_hemi_img,surf_data_hemi)
            
            #Export preproceced data
            out_flt_file_hemi = "{}/{}_{}.func.gii".format(pp_data_func_dir,func_fn_hemi.split('/')[-1][:-9],high_pass_type) 
            nb.save(flt_im_hemi, out_flt_file_hemi)
            
        for task in tasks:
            print(task)
     
            # preproc tasks runs
            preproc_files_hemi = glob.glob("{}/*_task-{}_*_hemi-{}_space-fsnative_bold_{}.func.gii".format(pp_data_func_dir,task,hemi, high_pass_type))

            
            
            # make correlation directory
            corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/fsnative".format(main_dir, project_dir, subject)
            os.makedirs(corr_dir, exist_ok=True)
    
            # output files 
            cor_file_hemi = "{}/{}_task-{}_hemi-{}_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,hemi,high_pass_type)

            # load preproc files to have meta and header
            preproc_img_hemi = nb.load(preproc_files_hemi[0])
            preproc_data_hemi = [x.data for x in preproc_img_hemi.darrays]
            preproc_data_hemi = np.vstack(preproc_data_hemi)
            
            preproc_img_header_hemi = preproc_img_hemi.header
            preproc_img_meta_hemi = preproc_img_hemi.meta
            
            # Correlations
            print('starting correlations')
            
            # compute the combination 
            combis_hemi = list(it.combinations(preproc_files_hemi, 2))
            
            #  load data and comute the correaltions
            cor_final_hemi = np.zeros((preproc_data_hemi.shape[1]))
            for combi_hemi in combis_hemi:
                task_cor_hemi = np.zeros((preproc_data_hemi.shape[1]))
                
                a_img_hemi = nb.load(combi_hemi[0])
                a_data_hemi = [x.data for x in a_img_hemi.darrays]
                a_data_hemi = np.vstack(a_data_hemi)
        
                b_img_hemi = nb.load(combi_hemi[1])
                b_data_hemi = [x.data for x in b_img_hemi.darrays]
                b_data_hemi = np.vstack(b_data_hemi)
                
                for vertice_hemi in range(a_data_hemi.shape[1]):
                    corr_hemi, _ = stats.pearsonr(a_data_hemi[:, vertice_hemi], b_data_hemi[:, vertice_hemi])
                    task_cor_hemi[vertice_hemi] = corr_hemi
            
                cor_final_hemi += task_cor_hemi / len(combis_hemi)
            
            print('starting exportation')
            
            # export correlations file
            
            cor_final_hemi = cor_final_hemi.astype(np.float32) 
            corr_img_hemi = nb.gifti.GiftiImage(header=preproc_img_header_hemi, meta=preproc_img_meta_hemi)
            corr_darray_hemi = nb.gifti.GiftiDataArray(data=cor_final_hemi.flatten())
            corr_img_hemi.add_gifti_data_array(corr_darray_hemi)
            nb.save(corr_img_hemi, cor_file_hemi)
        
  
    
            # Averaging 
            print('averaging')
            # make avg dir
            avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_avg/fsnative".format(main_dir, project_dir, subject)
            os.makedirs(avg_dir, exist_ok=True)
            
            # avg files 
            avg_file_hemi = "{}/{}_task-{}_hemi-{}_fmriprep_{}_avg_bold.func.gii".format(avg_dir, subject, task,hemi,high_pass_type)
 
            print("averaging...")
            data_avg_hemi = np.zeros(preproc_data_hemi.shape)
            
            for preproc_file_hemi in preproc_files_hemi:
                
                print('add: {}'.format(preproc_file_hemi))

                data_val_hemi = []
                # Load data
                preproc_img_hemi = nb.load(preproc_file_hemi)

                
                data_val_hemi = [x.data for x in preproc_img_hemi.darrays]
                data_val_hemi = np.vstack(data_val_hemi)
                
                # Averaging 
                data_avg_hemi += data_val_hemi/len(preproc_files_hemi)

            # export averaging data
            avg_img_hemi = make_giti_image(preproc_img_hemi,data_avg_hemi)
            nb.save(avg_img_hemi, avg_file_hemi)
                
                
 
            # Leave-one-out averages
            if len(preproc_files_hemi):
                combi_hemi = list(it.combinations(preproc_files_hemi, len(preproc_files_hemi)-1))
            
            for loo_num, avg_runs in enumerate(combi_hemi):
                loo_avg_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_loo_avg/fsnative".format(main_dir, project_dir, subject)
                os.makedirs(loo_avg_dir, exist_ok=True)
                
                print("loo_avg-{}".format(loo_num+1))
                
                # compute average between loo runs
                loo_avg_file_hemi = "{}/{}_task-{}_hemi-{}_fmriprep_{}_avg_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,hemi,high_pass_type, loo_num+1)
    
                # load data and make the loo_avg object
                preproc_img_hemi = nb.load(preproc_files_hemi[0])
                
                preproc_data_hemi = [x.data for x in preproc_img_hemi.darrays]
                preproc_data_hemi = np.vstack(preproc_data_hemi)
                
                data_loo_avg_hemi = np.zeros(preproc_data_hemi.shape)
            
                # compute leave on out averagin
                for avg_run in avg_runs:
                    print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
                    data_val_hemi = []
                    data_val_img_hemi = nb.load(avg_run)
                    data_val_hemi = [x.data for x in data_val_img_hemi.darrays]
                    data_val_hemi = np.vstack(data_val_hemi)
                    data_loo_avg_hemi += data_val_hemi/len(avg_runs)
                    
                
                
                # export leave one out file 
                loo_avg_img_hemi = make_giti_image(preproc_img_hemi,data_loo_avg_hemi)
                

                                                          
                nb.save(loo_avg_img_hemi, loo_avg_file_hemi)        
            
                # copy loo run (left one out run)
                for loo in preproc_files_hemi:
                    if loo not in avg_runs:
                        loo_file_hemi =  "{}/{}_task-{}_hemi-{}_fmriprep_{}_loo-{}_bold.func.gii".format(loo_avg_dir, subject,task,hemi,high_pass_type, loo_num+1)
                        print("loo: {}".format(loo))
                        os.system("{} {} {}".format(trans_cmd, loo, loo_file_hemi))

# Anatomy
print("getting anatomy...")
orig_dir_anat = "{}/{}/derivatives/fmriprep/fmriprep/{}/ses-01/anat/".format(main_dir, project_dir, subject)
pycortex_flat_dir = '{}/{}/derivatives/pp_data/cortex/db/{}/surfaces'.format(main_dir,project_dir,subject)
anat_files = glob.glob("{}/*.surf.gii".format(orig_dir_anat))

dest_dir_anat = "{}/{}/derivatives/pp_data/{}/anat".format(main_dir, project_dir, subject)
os.makedirs(dest_dir_anat, exist_ok=True)

for hemi in hemis : 
    if hemi == 'L' :
        flat_img_l = nb.load('{}/flat_lh.gii'.format(pycortex_flat_dir))
        flat_img_l.darrays[0].meta['AnatomicalStructurePrimary'] = 'CortexLeft'
        flat_img_l.darrays[0].meta['GeometricType']= 'Flat'
        nb.save(flat_img_l,'{}/{}_flat_lh.surf.gii'.format(dest_dir_anat,subject))
        
    elif hemi == 'R' :
        flat_img_r = nb.load('{}/flat_rh.gii'.format(pycortex_flat_dir))
        flat_img_r.darrays[0].meta['AnatomicalStructurePrimary'] = 'CortexRight'
        flat_img_r.darrays[0].meta['GeometricType']= 'Flat'
        nb.save(flat_img_r,'{}/{}_flat_rh.surf.gii'.format(dest_dir_anat,subject))

for orig_file in anat_files:
    file_name = os.path.basename(orig_file)
    dest_file = os.path.join(dest_dir_anat, file_name)
    shutil.copyfile(orig_file, dest_file)
    
            

    
    
