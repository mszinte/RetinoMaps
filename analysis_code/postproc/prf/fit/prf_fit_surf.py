"""
-----------------------------------------------------------------------------------------
prf_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
pRF fit code run by submit_fit.py
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: input filepath (timeseries)
sys.argv[2]: visual design filepath
sys.argv[3]: model output filepath
sys.argv[4]: timeseries prediction output path
sys.argv[5]: number of processors
-----------------------------------------------------------------------------------------
Output(s):
Nifti image files with fit parameters for a z slice
-----------------------------------------------------------------------------------------
To run :
>> cd to function directory
cd ~/projects/amblyo_prf/analysis_code/postproc/
>> python prf/fit/prf_fit.py [subject] [timeseries] [visual design] 
                     [fit] [prediction] [nb_procs]
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import sys, os
import numpy as np
import glob
import datetime
import json
import ipdb
deb = ipdb.set_trace

# MRI analysis imports
# from prfpy.rf import *
# from prfpy.timecourse import *
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter
import nibabel as nb

# Get inputs
start_time = datetime.datetime.now()
subject = sys.argv[1]
input_fn_HCP = sys.argv[2]
input_fn_fsnative = sys.argv[3]
input_vd = sys.argv[4]
fit_fn_HCP = sys.argv[5]
fit_fn_fsnative = sys.argv[5]
pred_fn_HCP = sys.argv[6]
pred_fn_fsnative = sys.argv[7]
nb_procs = int(sys.argv[8])

# Analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
screen_size_cm = analysis_info['screen_size_cm']
screen_distance_cm = analysis_info['screen_distance_cm']
TR = analysis_info['TR']
grid_nr = analysis_info['grid_nr']
max_ecc_size = analysis_info['max_ecc_size']

# Get task specific visual design matrix
vdm = np.load(input_vd)

# Load HCP 170k data 
data_img_HCP = nb.load(input_fn_HCP)
data_HCP = data_img_HCP.get_fdata()

data_var_HCP = np.var(data_HCP,axis=-1)
mask_HCP = data_var_HCP!=0.0    
num_vox_HCP = mask_HCP[...].sum()
data_to_analyse_HCP = data_HCP[mask_HCP]
data_where_HCP = np.where(data_var_HCP!=0.0)
data_indices_HCP = []

for x,y in zip(data_where_HCP[0],data_where_HCP[1]):
    data_indices_HCP.append((x,y))
fit_mat_HCP = np.zeros((data_HCP.shape[0],data_HCP.shape[1],6))
pred_mat_HCP = np.zeros(data_HCP.shape)

# determine model
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

gauss_model = Iso2DGaussianModel(stimulus=stimulus)
sizes = max_ecc_size * np.linspace(0.1,1,grid_nr)**2
eccs = max_ecc_size * np.linspace(0.25,1,grid_nr)**2
polars = np.linspace(0, 2*np.pi, grid_nr)

# grid fit
print("Grid fit")
gauss_fitter = Iso2DGaussianFitter(data=data_to_analyse_HCP, model=gauss_model, n_jobs=nb_procs)
gauss_fitter.grid_fit(ecc_grid=eccs, polar_grid=polars, size_grid=sizes)

# iterative fit
print("Iterative fit")
gauss_fitter.iterative_fit(rsq_threshold=0.0001, verbose=False)
fit_fit_HCP = gauss_fitter.iterative_search_params

# Re-arrange data
for est,vox in enumerate(data_indices_HCP):
    fit_mat_HCP[vox] = fit_fit_HCP[est]
    pred_mat_HCP[vox] = gauss_model.return_prediction(  mu_x=fit_fit_HCP[est][0], mu_y=fit_fit_HCP[est][1], size=fit_fit_HCP[est][2], 
                                                    beta=fit_fit_HCP[est][3], baseline=fit_fit_HCP[est][4])



fit_img_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=fit_mat_HCP, header=data_img_HCP.header)
nb.save(fit_img_HCP, fit_fn_HCP)

pred_img_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=pred_mat_HCP, header=data_img_HCP.header)
nb.save(pred_img_HCP, pred_fn_HCP)

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))


