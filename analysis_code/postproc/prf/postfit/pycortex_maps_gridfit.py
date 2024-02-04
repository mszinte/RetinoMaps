"""
-----------------------------------------------------------------------------------------
pycortex_maps_gridfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_gridfit.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_maps_gridfit.py ~/disks/meso_shared RetinoMaps sub-20 y
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import cortex
import importlib
import ipdb
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file,load_surface_pycortex

deb = ipdb.set_trace

#Define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg_in = sys.argv[4]

try:
    if save_svg_in == 'yes' or save_svg_in == 'y':
        save_svg = True
    elif save_svg_in == 'no' or save_svg_in == 'n':
        save_svg = False
    else:
        raise ValueError
except ValueError:
    sys.exit('Error: incorrect input (Yes, yes, y or No, no, n)')
       
# Maps settings
rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
    amp_idx, baseline_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8
      
cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

description_end = 'avg gridfit'
deriv_fn_label = 'avg-gridfit'

# plot scales
rsq_scale = [0, 1]
ecc_scale = [0, 7.5]
size_scale = [0, 7.5]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)
 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject,format_)
    fit_dir = "{}/fit".format(prf_dir)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    flatmaps_dir = '{}/pycortex/flatmaps_avg_gauss_gridfit'.format(prf_dir)
    datasets_dir = '{}/pycortex/datasets_avg_gauss_gridfit'.format(prf_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    if format_ == 'fsnative':
        deriv_avg_fn_L = '{}/{}_task-pRF_hemi-L_fmriprep_dct_avg_prf-deriv_gauss_gridfit.func.gii'.format(prf_deriv_dir, subject)
        deriv_avg_fn_R = '{}/{}_task-pRF_hemi-R_fmriprep_dct_avg_prf-deriv_gauss_gridfit.func.gii'.format(prf_deriv_dir, subject)
        deriv_mat = load_surface_pycortex(L_fn=deriv_avg_fn_L, R_fn=deriv_avg_fn_R)
        
    elif format_ == '170k':
        deriv_avg_fn = '{}/{}_task-pRF_fmriprep_dct_avg_prf-deriv_gauss_gridfit.dtseries.nii'.format(prf_deriv_dir, subject)
        deriv_mat = load_surface_pycortex(brain_fn=deriv_avg_fn)
        save_svg = False
    
    print('Creating flatmaps...')
    
    maps_names = []
    
    # threshold data
    deriv_mat_th = deriv_mat
    amp_down =  deriv_mat_th[amp_idx,...] > 0
    rsqr_th_down = deriv_mat_th[rsq_idx,...] >= analysis_info['rsqr_th'][0]
    rsqr_th_up = deriv_mat_th[rsq_idx,...] <= analysis_info['rsqr_th'][1]
    size_th_down = deriv_mat_th[size_idx,...] >= analysis_info['size_th'][0]
    size_th_up = deriv_mat_th[size_idx,...] <= analysis_info['size_th'][1]
    ecc_th_down = deriv_mat_th[ecc_idx,...] >= analysis_info['ecc_th'][0]
    ecc_th_up = deriv_mat_th[ecc_idx,...] <= analysis_info['ecc_th'][1]
    all_th = np.array((amp_down,rsqr_th_down,rsqr_th_up,size_th_down,size_th_up,ecc_th_down,ecc_th_up)) 
    deriv_mat[rsq_idx,np.logical_and.reduce(all_th)==False]=0
    
    # r-square
    rsq_data = deriv_mat[rsq_idx,...]
    alpha_range = analysis_info["alpha_range"]
    alpha = (rsq_data - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha[alpha>1]=1
    param_rsq = {'data': rsq_data, 'cmap': cmap_uni, 'alpha': rsq_data, 
                 'vmin': rsq_scale[0], 'vmax': rsq_scale[1], 'cbar': 'discrete', 
                 'cortex_type': 'VertexRGB','description': 'pRF rsquare',
                 'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg,
                 'cbar_label': 'pRF R2', 'with_labels': True}
    maps_names.append('rsq')
    
    # polar angle
    pol_comp_num = deriv_mat[polar_real_idx,...] + 1j * deriv_mat[polar_imag_idx,...]
    polar_ang = np.angle(pol_comp_num)
    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)
    ang_norm = np.fmod(ang_norm + col_offset,1)
    param_polar = {'data': ang_norm, 'cmap': cmap_polar, 'alpha': alpha, 
                   'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps, 'cortex_type': 'VertexRGB',
                   'cbar': 'polar', 'col_offset': col_offset, 
                   'description': 'pRF polar:{:3.0f} steps{}'.format(cmap_steps, description_end), 
                   'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg, 
                   'with_labels': True}
    exec('param_polar_{cmap_steps} = param_polar'.format(cmap_steps = int(cmap_steps)))
    exec('maps_names.append("polar_{cmap_steps}")'.format(cmap_steps = int(cmap_steps)))
    
    # eccentricity
    ecc_data = deriv_mat[ecc_idx,...]
    param_ecc = {'data': ecc_data, 'cmap': cmap_ecc_size, 'alpha': alpha,
                 'vmin': ecc_scale[0], 'vmax': ecc_scale[1], 'cbar': 'ecc', 'cortex_type': 'VertexRGB',
                 'description': 'pRF eccentricity{}'.format(description_end), 'curv_brightness': 1,
                 'curv_contrast': 0.1, 'add_roi': save_svg, 'with_labels': True}
    
    maps_names.append('ecc')
    
    # size
    size_data = deriv_mat[size_idx,...]
    param_size = {'data': size_data, 'cmap': cmap_ecc_size, 'alpha': alpha, 
                  'vmin': size_scale[0], 'vmax': size_scale[1], 'cbar': 'discrete', 
                  'cortex_type': 'VertexRGB', 'description': 'pRF size{}'.format(description_end), 
                  'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': False, 'cbar_label': 'pRF size',
                  'with_labels': True}
    maps_names.append('size')
    
    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = 'pRF_{}'.format(maps_name)
        roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
        exec("plt.savefig('{}/{}_task-pRF_{}_{}.pdf')".format(flatmaps_dir, subject, maps_name, deriv_fn_label))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_task-pRF_{}.hdf".format(datasets_dir, subject, deriv_fn_label)
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)
    
    
