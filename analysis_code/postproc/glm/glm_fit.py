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
python glm_fit.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
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
import nibabel as nb
import matplotlib.pyplot as plt

# nilearn import
from nilearn.plotting import plot_design_matrix 
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
tasks = analysis_info['task_glm']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
confounds_list = analysis_info['glm_confounds']
func_session = analysis_info['func_session'][0]

for format_, extension in zip(formats, extensions):
    # make folders
    glm_dir = '/{}/{}/derivatives/pp_data/{}/{}/glm/glm_fit'.format(main_dir, 
                                                                    project_dir, 
                                                                    subject, 
                                                                    format_)
    os.makedirs(glm_dir, exist_ok=True)

    for task in tasks :             
        # prepoc files name
        preproc_fns = glob.glob('{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct_loo_avg/*task-{}*dct_avg_loo*.{}'.format(main_dir, project_dir, subject, format_, task, extension))

        for preproc_fn in preproc_fns :
            match = re.search(r'_loo-(\d+)_', preproc_fn)
            loo_num = 'loo-{}'.format(match.group(1))
        
        
        

        
            
            # find the events and confounds files 
            event_dir = '{}/{}/{}/{}/func/'.format(main_dir, 
                                                   project_dir, 
                                                   subject, 
                                                   func_session)
            
            con_dir = '{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func'.format(main_dir, 
                                                                       project_dir, 
                                                                       subject, 
                                                                       func_session)
            
                
        

            # # Find the event files 
            event_file = glob.glob("{}/{}_{}_task-{}_run-01_events.tsv".format(event_dir, 
                                                                              subject, 
                                                                              func_session, 
                                                                              task))

        
            # make the designe matrixe  
            events = eventsMatrix(design_file=event_file[0], task=task, tr=TR)
            
            frame_times = np.arange(TRs) * TR
            design_matrix = make_first_level_design_matrix(frame_times,
                                                       events=events,
                                                       hrf_model='spm',
                                                       drift_model=None)
            
            design_matrix = design_matrix.drop(columns='Fix')
            #  Save the designe matrix 
            dm_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/designe_matrix'.format(main_dir, 
                                                                                 project_dir, 
                                                                                 subject, 
                                                                                 format_)
            os.makedirs(dm_dir, exist_ok=True)
            design_matrix.to_csv('{}/{}_task-{}_designe_matrix.tsv'.format(dm_dir, subject, task), sep="\t", na_rep='NaN', index=False)
            plt.figure()
            plot_design_matrix(design_matrix)
            plt.savefig('{}/{}_task-{}_designe_matrix.pdf'.format(dm_dir, subject, task))
            plt.close()

            # make glm output filenames
            glm_pred_fn = preproc_fn.split('/')[-1].replace('bold', 'glm-pred') 
            glm_fit_fn = preproc_fn.split('/')[-1].replace('bold', 'glm-fit') 

            # Load data
            preproc_img, preproc_data = load_surface(fn=preproc_fn)

            # Run the glm
            labels, estimates = run_glm(preproc_data, design_matrix.values, noise_model="ols")
            
            # extract glm predictions and r2       
            glm_pred, glm_r2 = extract_predictions_r2 (labels=labels,
                                                        estimate=estimates,
                                                        source_data=preproc_data)
        
            
            # export pred
            pred_img = make_surface_image(data=glm_pred, 
                                          source_img=preproc_img)

            nb.save(pred_img,'{}/{}'.format(glm_dir, glm_pred_fn)) 
             
            print('{} is done'.format(preproc_fn))


        
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))

        
# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))

