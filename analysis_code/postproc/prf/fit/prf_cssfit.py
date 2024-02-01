"""
-----------------------------------------------------------------------------------------
prf_cssfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing css fit
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
python prf_cssfit.py [main directory] [project name] [subject name] 
[inout file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_cssfit.py /scratch/mszinte/data RetinoMaps sub-02 /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-02/fsnative/func/fmriprep_dct_loo_avg/sub-02_task-pRF_hemi-L_fmriprep_dct_avg_loo-4_bold.func.gii 32  
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
import importlib
import numpy as np
deb = ipdb.set_trace


# MRI analysis imports
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel, CSS_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, CSS_Iso2DGaussianFitter
import nibabel as nb
import cortex


# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import load_surface ,make_surface_image
from pycortex_utils import set_pycortex_config_file, get_roi_verts_hemi
from maths_utils import r2_score_surf

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
css_params_num = 9


# Analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
vdm_width = analysis_info['vdm_size_pix'][0] 
vdm_height = analysis_info['vdm_size_pix'][1]
TR = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
rois = analysis_info['rois']
rsq_iterative_th = analysis_info['rsq_iterative_th']
css_grid_nr = analysis_info['css_grid_nr']

# Define directories and files names (fn)
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/prf/fit".format(
        main_dir,project_dir,subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/prf/fit".format(
        main_dir,project_dir,subject)
    os.makedirs(prf_fit_dir, exist_ok=True)


fit_fn_css = input_fn.split('/')[-1]
fit_fn_css = fit_fn_css.replace('bold', 'prf-fit_css')

pred_fn_css = input_fn.split('/')[-1]
pred_fn_css = pred_fn_css.replace('bold', 'prf-pred_css')

vdm_fn = "{}/{}/derivatives/vdm/vdm_pRF_{}_{}.npy".format(main_dir, project_dir, vdm_width, vdm_height)
pycortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)

# Set pycortex db and colormaps
set_pycortex_config_file(pycortex_dir)
importlib.reload(cortex)

# Get task specific visual design matrix
vdm = np.load(vdm_fn)

# defind model parameter grid range
sizes = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
polars = np.linspace(0, 2*np.pi, gauss_grid_nr)
exponent_css_grid = np.linspace(0, 4, css_grid_nr)


# load data
img, data, data_roi, roi_idx = get_roi_verts_hemi(fn=input_fn,subject=subject,rois=rois)



print('roi extraction done')





# determine visual design
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

# Gauss fit
# -------

# define gauss model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data_roi.T, 
                                   model=gauss_model, 
                                   n_jobs=n_jobs)



gauss_fitter.grid_fit(ecc_grid=eccs, 
                      polar_grid=polars, 
                      size_grid=sizes, 
                      verbose=verbose, 
                      n_batches=n_batches)


# iterative fit gauss model
gauss_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, verbose=verbose)
gauss_fit = gauss_fitter.iterative_search_params

# CSS fit
# -------

# define CSS model
css_model = CSS_Iso2DGaussianModel(stimulus=stimulus)

# grid fit CSS model
css_fitter = CSS_Iso2DGaussianFitter(data=data_roi.T, 
                                     model=css_model, 
                                     n_jobs=n_jobs,
                                     use_previous_gaussian_fitter_hrf=False,
                                     previous_gaussian_fitter=gauss_fitter
                                    )



# run css grid fit
css_fitter.grid_fit(exponent_grid=exponent_css_grid,
                    verbose=verbose,
                    n_batches=n_batches)

# run iterative fit
css_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, 
                         verbose=verbose, 
                         xtol=1e-4, 
                         ftol=1e-4)

css_fit = css_fitter.iterative_search_params

# rearange result of CSS model 
css_fit_mat = np.full((data.shape[1], css_params_num), np.nan, dtype=float)
css_pred_mat = np.full_like(data, np.nan, dtype=float) 


for est,vert in enumerate(roi_idx):

    css_fit_mat[vert] = css_fit[est]
    css_pred_mat[:,vert] = css_model.return_prediction(mu_x=css_fit[est][0],
                                                      mu_y=css_fit[est][1], 
                                                      size=css_fit[est][2], 
                                                      beta=css_fit[est][3], 
                                                      baseline=css_fit[est][4],
                                                      n=css_fit[est][5],
                                                      hrf_1=css_fit[est][6],
                                                      hrf_2=css_fit[est][7])
        


# compute loo r2
loo_bold_fn = input_fn.replace('dct_avg_loo','dct_loo')
loo_img, loo_bold = load_surface(fn=loo_bold_fn)

loo_r2 = r2_score_surf(bold_signal=loo_bold, model_prediction=css_pred_mat)

# add loo r2 css_fit_mat
css_fit_mat = np.column_stack((css_fit_mat, loo_r2))



# export data
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
              'n', 'hrf_1','hrf_2', 'r_squared','loo_r2']



# export fit
img_css_fit_mat = make_surface_image(data=css_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_css_fit_mat,'{}/{}'.format(prf_fit_dir, fit_fn_css)) 

# export pred
img_css_pred_mat = make_surface_image(data=css_pred_mat, source_img=img)
nb.save(img_css_pred_mat,'{}/{}'.format(prf_fit_dir, pred_fn_css)) 


# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))










