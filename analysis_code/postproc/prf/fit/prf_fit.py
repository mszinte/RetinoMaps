"""
-----------------------------------------------------------------------------------------
prf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
High-pass filter, z-score, run correlations, average, loov average and pick anat files
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit
2. run python command
python preproc_end_surf.py [main directory] [project name] [subject name] [input_fn] [n_jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_fit.py /scratch/mszinte/data RetinoMaps sub-02 /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-02/func/fmriprep_dct_avg/fsnative/sub-02_task-pRF_hemi-L_fmriprep_dct_avg_bold.func.gii 32  
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
from prfpy.model import Iso2DGaussianModel,Norm_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
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
n_jobs = sys.argv[5]
n_batches = n_jobs
verbose = True


# Analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
TR = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
dn_grid_nr = analysis_info['dn_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
rsq_threshold = analysis_info['fit_rsq_threshold']

# Define directories


prf_fit_dir = "{}/{}/derivatives/pp_data/{}/prf/fit".format(
    main_dir,project_dir,subject)
os.makedirs(prf_fit_dir, exist_ok=True)



fit_fn_DN = input_fn.split('/')[-1]
fit_fn_DN = fit_fn_DN.replace('bold', 'prf-fit')

pred_fn_DN = input_fn.split('/')[-1]
pred_fn_DN = fit_fn_DN.replace('bold', 'prf-pred')



vdm_fn = "{}/{}/derivatives/vdm/vdm.npy".format(main_dir, project_dir)

# Get task specific visual design matrix
vdm = np.load(vdm_fn)

# defind pRF parameter range
sizes = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.25,1,gauss_grid_nr)**2
polars = np.linspace(0, 2*np.pi, gauss_grid_nr)

# defind dn parameters
fixed_grid_baseline = 0
grid_bounds = [(0,1000),(0,1000)]
surround_size_grid = max_ecc_size * np.linspace(0.1, 1, dn_grid_nr)**2
surround_amplitude_grid = np.linspace(0, 10, dn_grid_nr)
surround_baseline_grid = np.linspace(0, 10, dn_grid_nr)
neural_baseline_grid = np.linspace(0, 10, dn_grid_nr)





gauss_bounds = [(-1.5*max_ecc_size, 1.5*max_ecc_size),  # x
                (-1.5*max_ecc_size, 1.5*max_ecc_size),  # y
                (0.1, 1.5*max_ecc_size),
                (0, 1000),  # prf amplitude
                (0, 1000),# prf size
                (0,10),(0,0)]  #hrf bounds


norm_bounds = [(-1.5*max_ecc_size, 1.5*max_ecc_size),  # x
                (-1.5*max_ecc_size, 1.5*max_ecc_size),  # y
                (0.1, 1.5*max_ecc_size),  # prf size
                (0, 1000),  # prf amplitude
                (0, 1000),  # bold baseline
                (0, 1000),  # surround amplitude
                (0.1, 3*max_ecc_size),  # surround size
                (0, 1000),  # neural baseline
                (1e-6, 1000),# surround baseline
                (0,10),(0,0)]  #hrf bounds








# load data
img, data = load_surface(fn=input_fn)
data = data[:,0:100] ## subsample


# determine gauss model
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)
gauss_model = Iso2DGaussianModel(stimulus=stimulus)



print('start grid gauss')
# grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data.T, model=gauss_model, n_jobs=n_jobs)
gauss_fitter.grid_fit(ecc_grid=eccs, 
                      polar_grid=polars, 
                      size_grid=sizes, 
                      verbose=verbose, 
                      n_batches=n_batches)


print('grid gauss ended')



# iterative fit Gauss model 
gauss_fitter.iterative_fit(rsq_threshold=rsq_threshold, 
                           verbose=verbose,
                          bounds = gauss_bounds,
                           xtol=1e-4,
                           ftol=1e-4)
gauss_fit = gauss_fitter.iterative_search_params

print('iterativ gauss ended')




# rearange result of Gauss model 
gauss_fit_mat = np.zeros((data.shape[1],8))
gauss_pred_mat = np.zeros_like(data) 
for est in range(len(data.T)):
    gauss_fit_mat[est] = gauss_fit[est]
    gauss_pred_mat[:,est] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                          mu_y=gauss_fit[est][1], 
                                                          size=gauss_fit[est][2], 
                                                          beta=gauss_fit[est][3], 
                                                          baseline=gauss_fit[est][4])

print('start grid dn')    
# determine DN model
dn_model = Norm_Iso2DGaussianModel(stimulus=stimulus)
dn_fitter = Norm_Iso2DGaussianFitter(data=data.T, 
                                     model=dn_model, 
                                     n_jobs=n_jobs,
                                     use_previous_gaussian_fitter_hrf=verbose,
                                     previous_gaussian_fitter=gauss_fitter)

# grid fit DN model  
dn_fitter.grid_fit(
    fixed_grid_baseline=fixed_grid_baseline,
    grid_bounds=grid_bounds,
    surround_amplitude_grid=surround_amplitude_grid,
    surround_size_grid=surround_size_grid,             
    surround_baseline_grid=surround_baseline_grid,
    neural_baseline_grid=neural_baseline_grid,
    n_batches=n_batches,
    rsq_threshold=rsq_threshold,
    verbose=verbose)

print('grid DN ended')



# iterative fit DN model 
dn_fitter.iterative_fit(rsq_threshold=rsq_threshold, 
                        verbose=verbose,
                        bounds = norm_bounds,
                        xtol=1e-4,
                        ftol=1e-4)
fit_fit_dn = dn_fitter.iterative_search_params

# rearange result of DN model 
dn_fit_mat = np.zeros((data.shape[1],12))
dn_pred_mat = np.zeros_like(data) 
for est in range(len(data.T)):
    dn_fit_mat[est] = fit_fit_dn[est]
    dn_pred_mat[:,est] = dn_model.return_prediction(mu_x=fit_fit_dn[est][0], 
                                                    mu_y=fit_fit_dn[est][1], 
                                                    prf_size=fit_fit_dn[est][2], 
                                                    prf_amplitude=fit_fit_dn[est][3], 
                                                    bold_baseline=fit_fit_dn[est][4],
                                                    srf_amplitude=fit_fit_dn[est][5],
                                                    srf_size=fit_fit_dn[est][6],
                                                    neural_baseline=fit_fit_dn[est][7],
                                                    surround_baseline=fit_fit_dn[est][8]
                                               )

print('iterativ DN ended')

#export data from DN model fit
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
              'srf_amplitude', 'srf_size','neural_baseline', 'surround_baseline',
              'hrf_1','hrf_2', 'r_squared']
              



# export fit
img_dn_fit_mat = make_surface_image(data=dn_fit_mat, source_img=img)
nb.save(img_dn_fit_mat,'{}/{}'.format(prf_fit_dir,fit_fn_DN)) 

# export pred
img_dn_pred_mat = make_surface_image(data=dn_pred_mat, source_img=img)
nb.save(img_dn_fit_mat,'{}/{}'.format(prf_fit_dir,pred_fn_DN)) 



for map_num, map_name in enumerate(maps_names):
    os.system('wb_command -set-map-names {}/{} -map {} {}'.format(prf_fit_dir,fit_fn_DN, map_num+1, map_name))


# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))




