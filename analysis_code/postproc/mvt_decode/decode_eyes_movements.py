"""
-----------------------------------------------------------------------------------------
decode_eyes_movements.py
-----------------------------------------------------------------------------------------
Goal of the script:
decode eyes movements or not
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/mvt_decode/
2. run python command
python preproc_end.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python decode_eyes_movements.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""
# General import 
import os
import nibabel as nb
import sys


import numpy as np



import json
import glob
import pandas as pd
from scipy.stats import loguniform 


import warnings
warnings.filterwarnings("ignore")

#sklearn import


from sklearn.svm import SVC
from sklearn.metrics import accuracy_score


from sklearn.model_selection import RandomizedSearchCV




main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
TR = analysis_info['TR']
high_pass_threshold = analysis_info['high_pass_threshold'] 
high_pass_type = analysis_info['high_pass_type'] 
sessions = analysis_info['session']
hemis = ['L','R']


task = 'SacLoc'
session = 'ses-02'
runs = ['01','02']
hemis = ['L','R']


# Datas end events directories
event_dir = '{}/{}/{}/{}/func/'.format(main_dir,project_dir,subject,session)
fsnative_bold_dir = '{}/{}/derivatives/pp_data/{}/func/fmriprep_dct/fsnative'.format(main_dir,project_dir,subject)
score_out_dir = '{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_decode_scores'.format(main_dir,project_dir,subject)
os.makedirs(score_out_dir, exist_ok=True)

# make an object with the events
# dim 1 = runs / dim 2 = TRs / dim 3 = eyemov_amplitude, eyemov_direction
stacked_events = [] 
for n_run, run in enumerate(runs) :
    event_file = glob.glob("{}/{}_{}_task-{}_run-{}_events.tsv".format(event_dir,subject,session,task,run))
    event_data = pd.read_table(event_file[0])
    eyemov_amplitude = np.array(event_data.eyemov_amplitude)
    eyemov_direction = np.array(event_data.eyemov_direction)
    amp_dir = np.column_stack((eyemov_amplitude, eyemov_direction))
    stacked_events.append(amp_dir)  
    
events_runs = np.stack(stacked_events, axis=0)

# Convert events in binary 0 = no eyes mvt 1 = eyes mvt
binar_events_runs = np.where(events_runs[:,:,1] == 33, 0, 1)



# Load BOLD data 
stacked_data_L = []
stacked_data_R = []

for n_run, run in enumerate(runs) :
    print(run)
    data_name_L = "{}/{}_{}_task-{}_run-{}_hemi-L_space-fsnative_bold_dct.func.gii".format(fsnative_bold_dir,subject,session,task,run)
    data_file_L = glob.glob(data_name_L)
    
    data_name_R = "{}/{}_{}_task-{}_run-{}_hemi-R_space-fsnative_bold_dct.func.gii".format(fsnative_bold_dir,subject,session,task,run)
    data_file_R = glob.glob(data_name_R)
    
    img_L = nb.load(data_file_L[0])
    data_L = [x.data for x in img_L.darrays]
    data_L = np.vstack(data_L) 
    header_L = img_L.header
    meta_L = img_L.meta
    print('hemi L done')

    
    img_R = nb.load(data_file_R[0])
    data_R = [x.data for x in img_R.darrays]
    data_R = np.vstack(data_R) 
    header_R = img_L.header
    meta_R = img_L.meta

    stacked_data_L.append(data_L)
    stacked_data_R.append(data_R)
    print('hemi R done')
    

data_runs_L = np.stack(stacked_data_L, axis=0) 
data_runs_R = np.stack(stacked_data_R, axis=0) 


param_dist = {
    'C': loguniform(1e-4, 1),
    'kernel' : ["linear", "poly", "rbf", "sigmoid"], 
    'gamma' : ['scale', 'auto']}
 

for hemi in hemis : 
    train_accuracy_verts = []  
    test_accuracy_verts = [] 
    
    out_file_name_train = '{}/{}_task-{}_hemi-{}_fmriprep_mcl_train_score_bold.func.gii'.format(score_out_dir,subject,task,hemi)
    
    out_file_name_test = '{}/{}_task-{}_hemi_{}_fmriprep_mcl_test_score_bold.func.gii'.format(score_out_dir,subject,task,hemi)
    
    
    if hemi == 'L':
        data_runs_hemi = data_runs_L
        header = header_L
        meta = meta_L
    elif hemi == 'R' :
        data_runs_hemi = data_runs_R
        header = header_R
        meta = meta_R
        
    for vert in range(data_runs_hemi.shape[2]):

        # split datas 
        X_train = data_runs_hemi[0, :, vert].reshape(-1, 1)
        X_test = data_runs_hemi[1, :, vert].reshape(-1, 1)
        y_train = binar_events_runs[0, :]
        y_test = binar_events_runs[1, :]

        # radom search to optimise 
        random_search = RandomizedSearchCV(SVC(), param_distributions=param_dist, n_iter=150, scoring='accuracy', random_state=42)
        random_search.fit(X_train, y_train)

        best_model = random_search.best_estimator_

        best_y_pred_train= best_model.predict(X_train)
        best_y_pred_test = best_model.predict(X_test)

        best_accuracy_score_train = accuracy_score(best_y_pred_train, y_train)
        best_accuracy_score_test = accuracy_score(best_y_pred_test, y_test)

        train_accuracy_verts.append(best_accuracy_score_train)
        test_accuracy_verts.append(best_accuracy_score_test)

    train_accuracy_verts = np.array(train_accuracy_verts)
    test_accuracy_verts = np.array(test_accuracy_verts)

    train_accuracy_verts_img = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in train_accuracy_verts:
        train_accuracy_verts_darray = nb.gifti.GiftiDataArray(data)
        train_accuracy_verts_img.add_gifti_data_array(train_accuracy_verts_darray)

    nb.save(train_accuracy_verts_img, out_file_name_train)

    test_accuracy_verts_img = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in train_accuracy_verts:
        test_accuracy_verts_darray = nb.gifti.GiftiDataArray(data)
        test_accuracy_verts_img.add_gifti_data_array(test_accuracy_verts_darray)

    nb.save(test_accuracy_verts_img, out_file_name_test)





