"""
-----------------------------------------------------------------------------------------
pycortex_maps_glm.py
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
>> cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/glm
2. run python command
>> python pycortex_maps_glm.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_maps_glm.py ~/disks/meso_shared RetinoMaps sub-02 n
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
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file,load_surface_pycortex

deb = ipdb.set_trace

#Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
tasks = analysis_info['task_glm']
alpha_range = analysis_info["alpha_range"]

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
z_map_idx, z_p_map_idx, fdr_idx, fdr_p_map_idx, r2_idx = 0,1,2,3,4
     
cmap = 'RdBu_r'
col_offset = 1.0/14.0
cmap_steps = 255



# plot scales
z_map_scale = [-7, 7]
p_map_scale = [0, 1]


# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)
 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    glm_dir = "{}/{}/derivatives/pp_data/{}/{}/glm".format(main_dir, project_dir, subject,format_)

    glm_deriv_dir = "{}/glm_derivatives".format(glm_dir)
    flatmaps_dir = '{}/pycortex/flatmaps_glm'.format(glm_dir)
    datasets_dir = '{}/pycortex/datasets_glm'.format(glm_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    for task in tasks : 
    
        if format_ == 'fsnative': 
            deriv_glm_fn_L = '{}/{}_task-{}_hemi-L_space-fsnative_dct_glm-fit_avg.func.gii'.format(glm_deriv_dir, subject, task)
            deriv_glm_fn_R = '{}/{}_task-{}_hemi-R_space-fsnative_dct_glm-fit_avg.func.gii'.format(glm_deriv_dir, subject, task)
            deriv_mat = load_surface_pycortex(L_fn=deriv_glm_fn_L, R_fn=deriv_glm_fn_R)
            
        elif format_ == '170k':
            deriv_avg_fn = '{}/{}_task-{}_space-fsLR_den-170k_dct_glm-fit_avg.dtseries.nii'.format(glm_deriv_dir, subject, task)
            deriv_mat = load_surface_pycortex(brain_fn=deriv_avg_fn)
            save_svg = False
        
        print('Creating flatmaps...')
        
        maps_names = []
        
   
        
        # p-map
        p_data = deriv_mat[z_p_map_idx,...]

        alpha = (p_data - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha[alpha>1]=1
        param_p_map = {'data': p_data, 'cmap': cmap ,
                      'vmin': p_map_scale[0], 'vmax': p_map_scale[1], 'cbar': 'discrete', 
                      'cortex_type': 'Vertex', 'description': 'p-map', 
                      'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg, 'cbar_label': 'p-value',
                      'with_labels': True}
        maps_names.append('p_map')
        
        # z-map
        z_data = deriv_mat[z_map_idx,...]
        param_z_map = {'data': z_data, 'cmap': cmap,
                      'vmin': z_map_scale[0], 'vmax': z_map_scale[1], 'cbar': 'discrete', 
                      'cortex_type': 'Vertex', 'description': 'z-map', 
                      'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg, 'cbar_label': 'z',
                      'with_labels': True}
        maps_names.append('z_map')
        
        
        # fdr_p-map
        fdr_p_data = deriv_mat[fdr_p_map_idx,...]

        alpha2 = (fdr_p_data - alpha_range[0])/(alpha_range[1]-alpha_range[0])
        alpha2[alpha2>1]=1
        param_fdr_p_map = {'data': fdr_p_data, 'cmap': cmap ,
                      'vmin': p_map_scale[0], 'vmax': p_map_scale[1], 'cbar': 'discrete', 
                      'cortex_type': 'Vertex', 'description': 'p-map', 
                      'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg, 'cbar_label': 'p-value',
                      'with_labels': True}
        maps_names.append('fdr_p_map')
        
        # fdr_z-map
        fdr_z_data = deriv_mat[fdr_idx,...]
        param_fdr_z_map = {'data': fdr_z_data, 'cmap': cmap ,
                      'vmin': z_map_scale[0], 'vmax': z_map_scale[1], 'cbar': 'discrete', 
                      'cortex_type': 'Vertex', 'description': 'z-map', 
                      'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg, 'cbar_label': 'z',
                      'with_labels': True}
        maps_names.append('fdr_z_map')
        


        
        # draw flatmaps
        volumes = {}
        for maps_name in maps_names:
        
            # create flatmap
            roi_name = 'glm_{}'.format(maps_name)
            roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
            print(roi_name)
            exec('param_{}.update(roi_param)'.format(maps_name))
            exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
            exec("plt.savefig('{}/{}_task-{}_{}.pdf')".format(flatmaps_dir, subject, task ,maps_name))
            plt.close()
        
            # save flatmap as dataset
            exec('vol_description = param_{}["description"]'.format(maps_name))
            exec('volume = volume_{}'.format(maps_name))
            volumes.update({vol_description:volume})
        
        # save dataset
        dataset_file = "{}/{}_task-{}.hdf".format(datasets_dir, task, subject)
        dataset = cortex.Dataset(data=volumes)
        dataset.save(dataset_file)
    
    
