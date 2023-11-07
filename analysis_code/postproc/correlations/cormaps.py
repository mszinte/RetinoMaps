# Imports
import os
import numpy as np
import nibabel as nb
import pandas as pd
import glob
import itertools as it
import ipdb
from scipy.stats.stats import pearsonr
debug = ipdb.set_trace

# Define folders, subjects, tasks
base_dir = '/scratch/mszinte/data/gaze_exp'
pybest_dir = "{}/derivatives/pybest".format(base_dir)

subjects = ['sub-001','sub-002']
task_set1 = ['GazeCW', 'GazeCCW']
task_set2 = ['GazeColumns', 'GazeLines']
task_sets = [task_set1,task_set2]
sessions = ['ses-04','ses-05']

# get files
for subject in subjects:
    print('Subject: {}'.format(subject))
    for task_set, session in zip(task_sets, sessions):
        print('Tasks: {} vs. {}'.format(task_set[0],task_set[1]))
        
        # Define files of task1
        task1_data_fns = glob.glob('{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_task-{task}_run-*_space-T1w_desc-preproc_bold.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[0]))
        task1_run_num = len(task1_data_fns)
        task1_runs = ['task-{}_run-{}'.format(task_set[0],run_num+1) for run_num in np.arange(0,task1_run_num)]
        task1_loo_avg_groups = list(it.combinations(task1_runs, task1_run_num-1))

        # Define files of task2
        task2_data_fns = glob.glob('{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_task-{task}_run-*_space-T1w_desc-preproc_bold.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[1]))
        task2_run_num = len(task2_data_fns)
        task2_runs = ['task-{}_run-{}'.format(task_set[1],run_num+1) for run_num in np.arange(0,task2_run_num)]
        task2_loo_avg_groups = list(it.combinations(task2_runs, task2_run_num-1))

        # get comparisons
        tasks_cor_groups = list(it.product(task1_loo_avg_groups,task2_loo_avg_groups))

        # Load brain masks
        task1_mask_fn = '{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_task-{task}_space-T1w_desc-preproc_mask.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[0])
        task1_mask = nb.load(task1_mask_fn).get_fdata()
        task2_mask_fn = '{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_task-{task}_space-T1w_desc-preproc_mask.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[1])
        task2_mask = nb.load(task2_mask_fn).get_fdata()
        task_mask = (task1_mask+task2_mask)/2.0

        # Create timeseries correlations with leave one out averaging
        for tasks_cor_group_num, tasks_cor_group in enumerate(tasks_cor_groups):
            print('Leave-one-out averaging combination {} in progress...'.format(tasks_cor_group_num+1))

            # Task1: load and average 
            for task1_cor_group_num, task1_cor_group in enumerate(tasks_cor_group[0]):
                task1_data_val = []
                task1_fn = '{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_{task_cor}_space-T1w_desc-preproc_bold.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[0], task_cor=task1_cor_group)
                print('Task 1: load + avg {}'.format(task1_fn))
                task1_data_val = nb.load(task1_fn).get_fdata()
                if task1_cor_group_num == 0: task1_data_avg = np.zeros(nb.load(task1_fn).shape)
                task1_data_avg += task1_data_val/len(task1_cor_group)

            # Task2: load and average 
            for task2_cor_group_num, task2_cor_group in enumerate(tasks_cor_group[1]):
                task2_data_val = []
                task2_fn = '{dir_fn}/{sub}/{ses}/preproc/{sub}_{ses}_{task_cor}_space-T1w_desc-preproc_bold.nii.gz'.format(
                                dir_fn=pybest_dir, sub=subject, ses=session, task=task_set[0], task_cor=task2_cor_group)
                print('Task 2: load + avg {}'.format(task2_fn))
                task2_data_val = nb.load(task2_fn).get_fdata()
                if task2_cor_group_num == 0: task2_data_avg = np.zeros(nb.load(task2_fn).shape)
                task2_data_avg += task2_data_val/len(task2_cor_group)

            
            # Mask data
            task1_data_avg_masked = task1_data_avg[task_mask==1,:]
            task2_data_avg_masked = task2_data_avg[task_mask==1,:]

            # Compute correlation
            print('Compute correlation of combination {} in progress...'.format(tasks_cor_group_num+1))
            tasks_cor = np.zeros(task1_data_avg_masked.shape[0])
            for voxnum in np.arange(0,task1_data_avg_masked.shape[0]):
                tasks_cor[voxnum] = pearsonr(task1_data_avg_masked[voxnum], task2_data_avg_masked[voxnum])[0]

            # Reshape dot product
            tasks_cor_reshape = np.zeros(task_mask.shape)
            tasks_cor_reshape[task_mask==1]=tasks_cor
            print("Combination {} correlation min: {:1.3f}, max: {:1.3f}\n".format(tasks_cor_group_num+1,np.min(tasks_cor_reshape),np.max(tasks_cor_reshape)))

            # Average corelations
            if tasks_cor_group_num == 0:tasks_cor_reshape_avg = np.zeros(tasks_cor_reshape.shape)
            tasks_cor_reshape_avg += tasks_cor_reshape/len(tasks_cor_groups)
        
        print("Average across {} combinations correlation min: {:1.3f}, max: {:1.3f}\n".format(len(tasks_cor_groups), np.min(tasks_cor_reshape_avg),np.max(tasks_cor_reshape_avg)))

        # Save correlation file
        cor_dir = "{}/derivatives/pp_data/{}/cor".format(base_dir,subject)
        try: os.makedirs(cor_dir)
        except: pass
        tasks_cor_fn = '{dir_fn}/{sub}_{task1}{task2}_timeseries_cor.nii.gz'.format(dir_fn=cor_dir, sub=subject, task1=task_set[0], task2=task_set[1])
        print('Saving correlation: {}'.format(tasks_cor_fn))
        tasks_cor_img = nb.Nifti1Image(dataobj=tasks_cor_reshape_avg, affine=nb.load(task1_mask_fn).affine, header=nb.load(task1_mask_fn).header)
        tasks_cor_img.to_filename(tasks_cor_fn)