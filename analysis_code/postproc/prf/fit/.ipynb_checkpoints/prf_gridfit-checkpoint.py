"""
-----------------------------------------------------------------------------------------
prf_gridfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing gaussian grid fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: number of jobs 
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit
2. run python command
python prf_gridfit.py [main directory] [project name] [subject name] 
[inout file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_gridfit.py /scratch/mszinte/data RetinoMaps sub-02 /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-02/func/fmriprep_dct_avg/fsnative/sub-02_task-pRF_hemi-L_fmriprep_dct_avg_bold.func.gii 32  
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import json
import ipdb
import datetime
import numpy as np
deb = ipdb.set_trace

# MRI analysis imports
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel 
from prfpy.fit import Iso2DGaussianFitter 
import nibabel as nb

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface

# Get inputs
start_time = datetime.datetime.now()


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
input_fn = sys.argv[4]
n_jobs = int(sys.argv[5])
n_batches = n_jobs
verbose = True
gauss_params_num = 8


# Analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
TR = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']


# Define directories
prf_fit_dir = "{}/{}/derivatives/pp_data/{}/prf/fit".format(
    main_dir,project_dir,subject)
os.makedirs(prf_fit_dir, exist_ok=True)

fit_fn_gauss_gridfit = input_fn.split('/')[-1]
fit_fn_gauss_gridfit = fit_fn_gauss_gridfit.replace('bold', 'prf-fit_gauss_gridfit')

pred_fn_gauss_gridfit = input_fn.split('/')[-1]
pred_fn_gauss_gridfit = pred_fn_gauss_gridfit.replace('bold', 'prf-pred_gauss_gridfit')

vdm_fn = "{}/{}/derivatives/vdm/vdm.npy".format(main_dir, project_dir)

# Get task specific visual design matrix
vdm = np.load(vdm_fn)

# defind model parameter grid range
sizes = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
polars = np.linspace(0, 2*np.pi, gauss_grid_nr)


# load data
img, data = load_surface(fn=input_fn)


# determine stimulus
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)


# determine gaussian model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data.T, 
                                   model=gauss_model, 
                                   n_jobs=n_jobs)

gauss_fitter.grid_fit(ecc_grid=eccs, 
                      polar_grid=polars, 
                      size_grid=sizes, 
                      verbose=verbose, 
                      n_batches=n_batches)



# rearange result of Gauss model 
gauss_fit = gauss_fitter.gridsearch_params
gauss_fit_mat = np.zeros((data.shape[1],gauss_params_num))
gauss_pred_mat = np.zeros_like(data) 



for est in range(len(data.T)):
    gauss_fit_mat[est] = gauss_fit[est]
    gauss_pred_mat[:,est] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                          mu_y=gauss_fit[est][1], 
                                                          size=gauss_fit[est][2], 
                                                          beta=gauss_fit[est][3], 
                                                          baseline=gauss_fit[est][4],
                                                          hrf_1=gauss_fit[est][5],
                                                          hrf_2=gauss_fit[est][6])



#export data from DN model fit
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline', 
              'hrf_1','hrf_2', 'r_squared']
              

# export fit
img_gauss_gridfit_fit_mat = make_surface_image(data=gauss_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_gauss_gridfit_fit_mat,'{}/{}'.format(prf_fit_dir, fit_fn_gauss_gridfit)) 

# export pred
img_gauss_gridfit_pred_mat = make_surface_image(data=gauss_pred_mat, source_img=img)
nb.save(img_gauss_gridfit_pred_mat,'{}/{}'.format(prf_fit_dir, pred_fn_gauss_gridfit)) 



# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))





