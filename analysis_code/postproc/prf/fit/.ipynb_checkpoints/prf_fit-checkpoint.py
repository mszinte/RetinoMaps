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
from prfpy.rf import *
from prfpy.timecourse import *
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter
import nibabel as nb

# Get inputs
start_time = datetime.datetime.now()
subject = sys.argv[1]
input_fn = sys.argv[2]
input_vd = sys.argv[3]
fit_fn = sys.argv[4]
pred_fn = sys.argv[5]
nb_procs = int(sys.argv[6])

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

# Load data
data_img = nb.load(input_fn)
data = data_img.get_fdata()

data_var = np.var(data,axis=-1)
mask = data_var!=0.0    
num_vox = mask[...].sum()
data_to_analyse = data[mask]
data_where = np.where(data_var!=0.0)
data_indices = []

for x,y,z in zip(data_where[0],data_where[1],data_where[2]):
    data_indices.append((x,y,z))
fit_mat = np.zeros((data.shape[0],data.shape[1],data.shape[2],6))
pred_mat = np.zeros(data.shape)

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
gauss_fitter = Iso2DGaussianFitter(data=data_to_analyse, model=gauss_model, n_jobs=nb_procs)
gauss_fitter.grid_fit(ecc_grid=eccs, polar_grid=polars, size_grid=sizes)

# iterative fit
print("Iterative fit")
gauss_fitter.iterative_fit(rsq_threshold=0.0001, verbose=False)
fit_fit = gauss_fitter.iterative_search_params

# Re-arrange data
for est,vox in enumerate(data_indices):
    fit_mat[vox] = fit_fit[est]
    pred_mat[vox] = gauss_model.return_prediction(  mu_x=fit_fit[est][0], mu_y=fit_fit[est][1], size=fit_fit[est][2], 
                                                    beta=fit_fit[est][3], baseline=fit_fit[est][4])

fit_img = nb.Nifti1Image(dataobj=fit_mat, affine=data_img.affine, header=data_img.header)
fit_img.to_filename(fit_fn)

pred_img = nb.Nifti1Image(dataobj=pred_mat, affine=data_img.affine, header=data_img.header)
pred_img.to_filename(pred_fn)

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))