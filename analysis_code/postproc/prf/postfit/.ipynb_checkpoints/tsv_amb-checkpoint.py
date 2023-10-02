#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 15:33:30 2023

@author: uriel
"""

## Imports
import os
import numpy as np
import pandas as pd
import nibabel as nb
import warnings
warnings.filterwarnings('ignore')


# Define parameters
n_subjects= 18 + 1
num = []
for t in range(1, n_subjects):
    num.append(str(t).zfill(2))
subjects = ['sub-{n}'.format(n=n) for n in num]
subjects_excluded = ['sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18'] # régler le problème avec sub-11
# subjects_plot = ['sub-001', 'sub-002', 'sub-003', 'sub-004',
#                  'sub-005', 'sub-006', 'sub-007', 'sub-008', 'group']

tasks = ['prf']
rois = ['V1','V2','V3','V3AB','LO','VO','hMT+','iIPS','sIPS']

# Define folders
base_dir = '/Users/uriel/disks/meso_shared/amblyo_prf'
bids_dir = "{}".format(base_dir)
pp_dir = "{}/derivatives/pp_data".format(base_dir)

# analysis settings
best_voxels_num = 250
type_analyses = ['','_best{}'.format(best_voxels_num)]

cortical_mask = 'cortical'
n_ecc_bins=10
verbose = False
TR = 1.8



### Compute TSV files for fullscreen atention R2 comparison


# Create TSV files
group_tsv_dir = '{}/{}/prf/tsv'.format(pp_dir, 'group')
try: os.makedirs(group_tsv_dir)
except: pass


for task in tasks:
    for subject_num, subject in enumerate(subjects):   
        if subject not in subjects_excluded[:]:
     
            print(subject_num)
            print(subject)
            # define folders
            fit_dir = '{}/{}/prf/fit'.format(pp_dir, subject)
            mask_dir = '{}/{}/masks'.format(pp_dir, subject)
            tsv_dir = '{}/{}/prf/tsv'.format(pp_dir, subject)
            try: os.makedirs(tsv_dir)
            except: pass
            
            # # load pRF threshold masks
            # th_mat = nb.load('{}/{}_task-{}_prf_threshold.nii.gz'.format(mask_dir,subject,task)).get_fdata()
    
            # load fit parameters x by threshold
            derives_mat = nb.load('{}/{}_task-prf_fmriprep_dct_bold_avg_prf-deriv.nii.gz'.format(fit_dir,subject)).get_fdata()
           
            df_rois = pd.DataFrame()
            # creat tsv
            for roi_num, roi in enumerate(rois):
                # load roi
                lh_mat = nb.load("{}/{}_{}_L.nii.gz".format(mask_dir, roi, cortical_mask)).get_fdata()
                rh_mat = nb.load("{}/{}_{}_R.nii.gz".format(mask_dir, roi, cortical_mask)).get_fdata()
                roi_mat = lh_mat + rh_mat
                roi_mat[roi_mat==0] = np.nan




                # select data by roi mask
                derives_roi_mat = derives_mat[roi_mat==True]
                
            
                # create dataframe
                df_roi = pd.DataFrame({'subject': [subject] * derives_roi_mat.shape[0],
                                       'roi': [roi] * derives_roi_mat.shape[0],
                                       'r2': derives_roi_mat[...,0],
                                       'ecc': derives_roi_mat[...,2],
                                       'sd': derives_roi_mat[...,5],
                                       'x': derives_roi_mat[...,8],
                                       'y': derives_roi_mat[...,9],
                                       'amplitude': derives_roi_mat[...,6],
                                       'baseline': derives_roi_mat[...,7]})
                
                df_rois = pd.concat([df_rois, df_roi], ignore_index=True)
                    
                
                
            
            # save dataframe
            df_fn = "{}/{}_task-{}_prf_threshold_par.tsv".format(tsv_dir,subject,task)
            print('saving {}'.format(df_fn))
            df_rois.to_csv(df_fn, sep="\t", na_rep='NaN',index=False)
            
            
            # across subject
            if subject_num == 0: df_group = df_rois
            else: df_group = pd.concat([df_group, df_rois])
            
        
        # save group data
        df_group_fn = "{}/group_task-{}_prf_threshold_par.tsv".format(group_tsv_dir,task)
        print('saving {}'.format(df_group_fn))
        df_group.to_csv(df_group_fn, sep="\t", na_rep='NaN')
        

