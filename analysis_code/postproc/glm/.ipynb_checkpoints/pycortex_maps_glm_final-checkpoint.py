"""
-----------------------------------------------------------------------------------------
pycortex_maps_glm_final.py
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
>> python pycortex_maps_glm_final.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_maps_glm_final.py ~/disks/meso_shared RetinoMaps sub-13 n
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
     
cmap = 'J4R'
col_offset = 1.0/14.0
cmap_steps = 4



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
    


    if format_ == 'fsnative': 
        deriv_glm_fn_L = '{}/{}_task-eyes-mvt_hemi-L_space-fsnative_dct_glm-significant_map.func.gii'.format(glm_deriv_dir, subject)
        deriv_glm_fn_R = '{}/{}_task-eyes-mvt_hemi-R_space-fsnative_dct_glm-significant_map.func.gii'.format(glm_deriv_dir, subject)
        final_mat = load_surface_pycortex(L_fn=deriv_glm_fn_L, R_fn=deriv_glm_fn_R)
        
    elif format_ == '170k':
        deriv_avg_fn = '{}/{}_task-eyes-mvt_space-fsLR_den-170k_dct_glm-significant_map.dtseries.nii'.format(glm_deriv_dir, subject)
        final_mat = load_surface_pycortex(brain_fn=deriv_avg_fn)
        save_svg = False
    
    print('Creating flatmaps...')
    
    maps_names = []
    
   
    
    # param_final = {'data': final_mat, 'cmap': cmap, '
    #                'vmin': 0, 'vmax': 3, 'cmap_steps': cmap_steps, 'cortex_type': 'Vertex',
    #                'cbar': 'polar', 'col_offset': col_offset, 
    #                'description': 'final map', 
    #                'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg, 
    #                'with_labels': True}
    
    # maps_names.append('final')
    
    final_data = final_mat[0,...]
    alpha_range = [0,0]
    alpha = (final_data - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha[alpha>1]=1
    param_final = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                 'vmin': 0, 'vmax': 3, 'cbar': 'discrete', 'cmap_steps': cmap_steps,
                 'cortex_type': 'VertexRGB','description': 'final map',
                 'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg,
                 'cbar_label': '', 'with_labels': True}
    maps_names.append('final')
    


    
    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = 'glm_{}'.format(maps_name)
        roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
        exec("plt.savefig('{}/{}_final_{}.pdf')".format(flatmaps_dir, subject, maps_name))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_final.hdf".format(datasets_dir, subject)
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)
    
    
