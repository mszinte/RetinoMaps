"""
-----------------------------------------------------------------------------------------
prf_fit_surf.py
-----------------------------------------------------------------------------------------
Goal of the script:
pRF fit adapted for surfaces, run by submit_fit.py
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
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit/
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
from prfpy.model import Iso2DGaussianModel,Norm_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
import nibabel as nb

# Get inputs
start_time = datetime.datetime.now()

subject = sys.argv[1]

if sys.argv[2].endswith('.nii'):
    input_fn_HCP = sys.argv[2]
elif sys.argv[2].endswith('.gii'):
    input_fn_fsnative = sys.argv[2]

input_vd = sys.argv[3]

if sys.argv[4].endswith('.nii'):
    fit_fn_HCP_gauss = sys.argv[4]
    fit_fn_HCP_DN = sys.argv[5]      
elif sys.argv[4].endswith('.gii'):
    fit_fn_fsnative_gauss = sys.argv[4]
    fit_fn_fsnative_DN = sys.argv[5] 

if sys.argv[6].endswith('.nii'):
    pred_fn_HCP_gauss = sys.argv[6]
    pred_fn_HCP_DN = sys.argv[7]
elif sys.argv[6].endswith('.gii'):
    pred_fn_fsnative_gauss = sys.argv[6]
    pred_fn_fsnative_DN = sys.argv[7]

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

# GIFTI
if 'input_fn_fsnative' in vars(): 

    # Load fsnative data 
    data_img_fsnative, data_fsnative = load_gifti_image(input_fn_fsnative)
    data_fsnative = data_fsnative.T

    fit_mat_gauss_fsnative = np.zeros((data_fsnative.shape[0],data_fsnative.shape[1],6))
    pred_mat_gauss_fsnative = np.zeros(data_fsnative.shape)
    
    fit_mat_DN_fsnative = np.zeros((data_fsnative.shape[0],data_fsnative.shape[1],6))
    pred_mat_DN_fsnative = np.zeros(data_fsnative.shape)
    
    # determine gauss model
    stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                             screen_distance_cm=screen_distance_cm,
                             design_matrix=vdm, 
                             TR=TR)
    
    gauss_model = Iso2DGaussianModel(stimulus=stimulus)
    sizes = max_ecc_size * np.linspace(0.1,1,grid_nr)**2
    eccs = max_ecc_size * np.linspace(0.25,1,grid_nr)**2
    polars = np.linspace(0, 2*np.pi, grid_nr)
    
    # grid fit
    print("Grid fit Gauss")
    gauss_fitter = Iso2DGaussianFitter(data=data_fsnative, model=gauss_model, n_jobs=nb_procs)
    gauss_fitter.grid_fit(ecc_grid=eccs, polar_grid=polars, size_grid=sizes, verbose=False, n_batches=8)
    
    # iterative fit
    print("Iterative fit Gauss")
    gauss_fitter.iterative_fit(rsq_threshold=0.0001, verbose=True)
    fit_fit_gauss_fsnative = gauss_fitter.iterative_search_params
    
    # determine DN model
    DN_model = Norm_Iso2DGaussianModel(stimulus=stimulus)
    
    DN_fitter = Norm_Iso2DGaussianFitter(data=data_fsnative, model=DN_model, n_jobs=nb_procs,
                                         use_previous_gaussian_fitter_hrf=True,
                                         previous_gaussian_fitter=gauss_fitter)
                                         
    
    
    # grid fit
    print("Grid fit DN")
    num = 7
    fixed_grid_baseline=0
    norm_grid_bounds = [(0,1000),(0,1000)] 
    
    DN_fitter.grid_fit(surround_amplitude_grid=np.linspace(0,10,num),
             surround_size_grid=np.linspace(1,10,num),             
             neural_baseline_grid=np.linspace(0,10,num),
             surround_baseline_grid=np.linspace(1,10,num),
             n_batches=8,
             rsq_threshold=0.0001,
             verbose = False,
             fixed_grid_baseline=fixed_grid_baseline,
             grid_bounds=norm_grid_bounds,
             hrf_1_grid=np.linspace(0,10,num),
             hrf_2_grid=np.linspace(0,0,1),
             # ecc_grid=eccs,
             # polar_grid=polars,
             # size_grid=sizes
             )
    
    
    
    
     
    # iterative fit
    print("Iterative fit Gauss")
    DN_fitter.iterative_fit(rsq_threshold=0.0001, verbose=False)
    fit_fit_DN_fsnative = DN_fitter.iterative_search_params
    


    
    print('start re arrange gauss')
    fit_mat_gauss_fsnative = fit_fit_gauss_fsnative
    # Re-arrange data
    for vert in range(len(fit_mat_gauss_fsnative)):
        pred_mat_gauss_fsnative[vert] = gauss_model.return_prediction(  mu_x=fit_fit_gauss_fsnative[vert][0], mu_y=fit_fit_gauss_fsnative[vert][1], size=fit_fit_gauss_fsnative[vert][2], 
                                                        beta=fit_fit_gauss_fsnative[vert][3], baseline=fit_fit_gauss_fsnative[vert][4])
    
    
    fit_mat_DN_fsnative = fit_fit_gauss_fsnative
    # Re-arrange data
    for vert in range(len(fit_mat_DN_fsnative)):

        pred_mat_DN_fsnative[vert] = gauss_model.return_prediction(  mu_x=fit_fit_DN_fsnative[vert][0], mu_y=fit_fit_DN_fsnative[vert][1], size=fit_fit_DN_fsnative[vert][2], 
                                                        beta=fit_fit_DN_fsnative[vert][3], baseline=fit_fit_DN_fsnative[vert][4])
    
    

    # export fsnative gauss fit 

    fit_img_gauss_fsnative = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in fit_mat_gauss_fsnative:
        darray = nb.gifti.GiftiDataArray(data)
        fit_img_gauss_fsnative.add_gifti_data_array(darray)
    nb.save(fit_img_gauss_fsnative, fit_fn_fsnative_gauss)
    
    
    # export fsnative gauss pred 
    pred_img_gauss_fsnative = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in pred_mat_gauss_fsnative:
        darray = nb.gifti.GiftiDataArray(data)
        pred_img_gauss_fsnative.add_gifti_data_array(darray)
    nb.save(pred_img_gauss_fsnative, pred_fn_fsnative_gauss)


    # export fsnative DN fit 
    fit_img_DN_fsnative = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in fit_mat_DN_fsnative:
        darray = nb.gifti.GiftiDataArray(data)
        fit_img_DN_fsnative.add_gifti_data_array(darray)
    nb.save(fit_img_DN_fsnative, fit_fn_fsnative_DN)
    
    # export fsnative DN pred
    pred_img_DN_fsnative = nb.gifti.GiftiImage(header = header, meta = meta)
    for data in pred_mat_DN_fsnative:
        darray = nb.gifti.GiftiDataArray(data)
        pred_img_DN_fsnative.add_gifti_data_array(darray)
    nb.save(pred_img_DN_fsnative, pred_fn_fsnative_DN)
    



# elif 'input_fn_HCP' in vars(): 


#     # Load HCP 170k data 
#     data_img_HCP = nb.load(input_fn_HCP)
#     data_HCP = data_img_HCP.get_fdata()
#     data_HCP = data_HCP[:,10:11 ] # juste for test
#     data_HCP = data_HCP.T # juste for test
    

        
        
        
#     fit_mat_gauss_HCP = np.zeros((data_HCP.shape[0],data_HCP.shape[1],6))
#     pred_mat_gauss_HCP = np.zeros(data_HCP.shape)
    
#     fit_mat_DN_HCP = np.zeros((data_HCP.shape[0],data_HCP.shape[1],6))
#     pred_mat_DN_HCP = np.zeros(data_HCP.shape)
    
#     # determine gauss model
#     stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
#                              screen_distance_cm=screen_distance_cm,
#                              design_matrix=vdm, 
#                              TR=TR)
    
#     gauss_model = Iso2DGaussianModel(stimulus=stimulus)
#     sizes = max_ecc_size * np.linspace(0.1,1,grid_nr)**2
#     eccs = max_ecc_size * np.linspace(0.25,1,grid_nr)**2
#     polars = np.linspace(0, 2*np.pi, grid_nr)
    
#     # grid fit
#     print("Grid fit Gauss")
#     gauss_fitter = Iso2DGaussianFitter(data=data_HCP, model=gauss_model, n_jobs=nb_procs)
#     gauss_fitter.grid_fit(ecc_grid=eccs, polar_grid=polars, size_grid=sizes, verbose=True)
    
 
#     # iterative fit
#     print("Iterative fit Gauss")
#     gauss_fitter.iterative_fit(rsq_threshold=0.0001, verbose=True)
#     fit_fit_gauss_HCP = gauss_fitter.iterative_search_params
    
    
    
#     # determine DN model
#     DN_model = Norm_Iso2DGaussianModel(stimulus=stimulus)
    
#     DN_fitter = Norm_Iso2DGaussianFitter(data=data_HCP, model=DN_model, n_jobs=nb_procs,
#                                          use_previous_gaussian_fitter_hrf=True,
#                                          previous_gaussian_fitter=gauss_fitter)
                                         
    
    
#     # grid fit
#     print("Grid fit DN")
#     num = 7
#     fixed_grid_baseline=0
#     norm_grid_bounds = [(0,1000),(0,1000)] 
#     DN_fitter.grid_fit(surround_amplitude_grid=np.linspace(0,10,num),
#              surround_size_grid=np.linspace(1,10,num),
#              neural_baseline_grid=np.linspace(0,10,num),
#              surround_baseline_grid=np.linspace(1,10,num),
#              fixed_grid_baseline=fixed_grid_baseline,
#              grid_bounds=norm_grid_bounds,
#              hrf_1_grid=np.linspace(0,10,num),
#              hrf_2_grid=np.linspace(0,0,1),
#              # ecc_grid=eccs,
#              # polar_grid=polars,
#              # size_grid=sizes
#              )
    
    
    
    
     
#     # iterative fit
#     print("Iterative fit Gauss")
#     DN_fitter.iterative_fit(rsq_threshold=0.0001, verbose=False)
#     fit_fit_DN_HCP = DN_fitter.iterative_search_params
    


    
#     print('start re arrange gauss')
#     fit_mat_gauss_HCP = fit_fit_gauss_HCP
#     # Re-arrange data
#     for vert in range(len(fit_mat_gauss_HCP)):
#         pred_mat_gauss_HCP[vert] = gauss_model.return_prediction(  mu_x=fit_fit_gauss_HCP[vert][0], mu_y=fit_fit_gauss_HCP[vert][1], size=fit_fit_gauss_HCP[vert][2], 
#                                                         beta=fit_fit_gauss_HCP[vert][3], baseline=fit_fit_gauss_HCP[vert][4])
    
    
#     fit_mat_DN_HCP = fit_fit_gauss_HCP
#     # Re-arrange data
#     for vert in range(len(fit_mat_DN_HCP)):

#         pred_mat_DN_HCP[vert] = gauss_model.return_prediction(  mu_x=fit_fit_DN_HCP[vert][0], mu_y=fit_fit_DN_HCP[vert][1], size=fit_fit_DN_HCP[vert][2], 
#                                                         beta=fit_fit_DN_HCP[vert][3], baseline=fit_fit_DN_HCP[vert][4])
    
    
#     print('export started',fit_fn_HCP)
#     # export HCP DN
#     fit_img_DN_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=fit_mat_DN_HCP, header=data_img_HCP.header)
#     nb.save(fit_img_DN_HCP, fit_fn_HCP)
    
#     print('export started',pred_fn_HCP)
#     pred_img_DN_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=pred_mat_DN_HCP, header=data_img_HCP.header)
#     nb.save(pred_img_DN_HCP, pred_fn_HCP)
    
#     print('export started',fit_fn_HCP)
#     # export HCP gauss
#     fit_img_gauss_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=fit_mat_gauss_HCP, header=data_img_HCP.header)
#     nb.save(fit_img_gauss_HCP, fit_fn_HCP)
    
    
#     print('export started',pred_fn_HCP)
#     pred_img_gauss_HCP = nb.cifti2.cifti2.Cifti2Image(dataobj=pred_mat_gauss_HCP, header=data_img_HCP.header)
#     nb.save(pred_img_gauss_HCP, pred_fn_HCP)
    
    

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))

