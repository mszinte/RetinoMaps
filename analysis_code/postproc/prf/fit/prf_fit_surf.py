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

from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel,Norm_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
import nibabel as nb


sys.path.append("{}/../../../utils".format(os.getcwd()))
from gifti_utils import make_gifti_image , load_gifti_image

# Get inputs
start_time = datetime.datetime.now()

subject = sys.argv[1]
input_vd = sys.argv[2]

if sys.argv[3].endswith('.nii'):
    input_fn_HCP = sys.argv[3]
    fit_fn_HCP_DN = sys.argv[4]
    pred_fn_HCP_DN = sys.argv[5]
    
    
elif sys.argv[3].endswith('.gii'):
    input_fn_fsnative = sys.argv[3]
    fit_fn_fsnative_DN = sys.argv[4] 
    pred_fn_fsnative_DN = sys.argv[5]
    
    
n_jobs = 8
n_batches = 8
rsq_threshold = 0.0001

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

# defind pRF parameter range
sizes = max_ecc_size * np.linspace(0.1,1,grid_nr)**2
eccs = max_ecc_size * np.linspace(0.25,1,grid_nr)**2
polars = np.linspace(0, 2*np.pi, grid_nr)


# defind dn parameters
fixed_grid_baseline = 0
grid_bounds = [(0,1000),(0,1000)]
surround_size_grid = sizes
surround_amplitude_grid = np.linspace(0, 10, grid_nr)
surround_baseline_grid = np.linspace(0, 10, grid_nr)
neural_baseline_grid = np.linspace(0, 10, grid_nr)

# GIFTI
if 'input_fn_fsnative' in vars(): 
    img_fsnative, data_fsnative = load_gifti_image(input_fn_fsnative)
    

    # determine gauss model
    stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1], 
                             screen_distance_cm=screen_distance_cm,
                             design_matrix=vdm, 
                             TR=TR)
    
    gauss_model = Iso2DGaussianModel(stimulus=stimulus)
    
    
    # grid fit gauss model
    print('gauss grid')
    gauss_fitter = Iso2DGaussianFitter(data=data_fsnative.T, model=gauss_model, n_jobs=n_jobs)
    gauss_fitter.grid_fit(ecc_grid=eccs, 
                          polar_grid=polars, 
                          size_grid=sizes, 
                          verbose=False, 
                          n_batches=n_batches)
    
    # iterative fit Gauss model 
    print('gauss iterative')
    gauss_fitter.iterative_fit(rsq_threshold=rsq_threshold, verbose=False)
    gauss_fit = gauss_fitter.iterative_search_params
    
    
    # rearange result of Gauss model 
    gauss_fit_mat = np.zeros((data_fsnative.shape[1],8))
    gauss_pred_mat = np.zeros_like(data_fsnative) 
    for est in range(len(data_fsnative.T)):
        gauss_fit_mat[est] = gauss_fit[est]
        gauss_pred_mat[:,est] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                              mu_y=gauss_fit[est][1], 
                                                              size=gauss_fit[est][2], 
                                                              beta=gauss_fit[est][3], 
                                                              baseline=gauss_fit[est][4])
    # determine DN model
    dn_model = Norm_Iso2DGaussianModel(stimulus=stimulus)
    dn_fitter = Norm_Iso2DGaussianFitter(data=data_fsnative.T, 
                                         model=dn_model, 
                                         n_jobs=n_jobs,
                                         use_previous_gaussian_fitter_hrf=True,
                                         previous_gaussian_fitter=gauss_fitter)
    

    # grid fit DN model  
    print('DN grid')
    dn_fitter.grid_fit(
        fixed_grid_baseline=fixed_grid_baseline,
        grid_bounds=grid_bounds,
        surround_amplitude_grid=surround_amplitude_grid,
        surround_size_grid=surround_size_grid,             
        surround_baseline_grid=surround_baseline_grid,
        neural_baseline_grid=neural_baseline_grid,
        n_batches=n_batches,
        rsq_threshold=rsq_threshold,
        verbose=False)
    
    # iterative fit DN model 
    print('DN grid')
    dn_fitter.iterative_fit(rsq_threshold=rsq_threshold, verbose=False)
    fit_fit_dn = dn_fitter.iterative_search_params
    
    
    # rearange result of DN model 
    dn_fit_mat = np.zeros((data_fsnative.shape[1],12))
    dn_pred_mat = np.zeros_like(data_fsnative) 
    for est in range(len(data_fsnative.T)):
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
        
    
    #export data from DN model fit
    maps_names = ['mu_x','mu_y','prf_size','prf_amplitude','bold_baseline',
                  'srf_amplitude','srf_size','surround_baseline ','hrf_1',
                  'hrf_2','r_squared']
    
    
    # export fit
    img_dn_fit_mat = make_gifti_image(img_fsnative,dn_fit_mat)
    nb.save(img_dn_fit_mat,fit_fn_fsnative_DN ) 
    
    # export pred
    img_dn_pred_mat = make_gifti_image(img_fsnative,dn_pred_mat)
    nb.save(img_dn_pred_mat,pred_fn_fsnative_DN ) 
    
    for map_num, mape_name in enumerate(maps_names):
        os.system('wb_command -set-map-names {} -map {} {}'.format(fit_fn_fsnative_DN,map_num+1,mape_name))





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

