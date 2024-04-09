"""
-----------------------------------------------------------------------------------------
pycortex_maps_stats_final.py
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
>> cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/stats
2. run python command
>> python pycortex_maps_glm_final.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_maps_stats_final.py ~/disks/meso_shared RetinoMaps sub-01 n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

#  Debug import 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import cortex
import importlib
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file,load_surface_pycortex

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
all_idx, pursuit_idx, saccade_idx, pursuit_and_saccade_idx, vision_idx, \
    vision_and_pursuit_idx, vision_and_saccade_idx, \
        vision_and_pursuit_and_saccade_idx = 0,1,2,3,4,5,6,7
     
cmap = 'stats_colors'
col_offset = 1.0/14.0
cmap_steps = 255


# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

if subject == 'sub-170k':
    formats = ['170k'] 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    stats_dir = "{}/{}/derivatives/pp_data/{}/{}/final_stats".format(main_dir, project_dir, subject,format_)
    stats_results_dir = "{}/results".format(stats_dir)
    
    flatmaps_dir = '{}/pycortex/flatmaps_stats'.format(stats_dir)
    datasets_dir = '{}/pycortex/datasets_stats'.format(stats_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    


    if format_ == 'fsnative': 
        deriv_stats_fn_L = '{}/{}_hemi-L_final-stats.func.gii'.format(stats_results_dir, subject)
        deriv_stats_fn_R = '{}/{}_hemi-R_final-stats.func.gii'.format(stats_results_dir, subject)
        
        results = load_surface_pycortex(L_fn=deriv_stats_fn_L, R_fn=deriv_stats_fn_R)
        final_mat = results['data_concat']
        
    elif format_ == '170k':
        stats_avg_fn = '{}/{}_final-stats.dtseries.nii'.format(stats_results_dir, subject)
        results = load_surface_pycortex(brain_fn=stats_avg_fn)
        final_mat = results['data_concat']
        
        if subject == 'sub-170k':
            save_svg = save_svg
        else: 
            save_svg = False
    
    print('Creating flatmaps...')
    
    maps_names = []
    
    #  Creat the all flatmap
    final_data = final_mat[all_idx,...]
    param_all = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                 'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                 'cortex_type': 'VertexRGB','description': 'final map',
                 'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                 'cbar_label': '', 'with_labels': True}
    maps_names.append('all')
    
    #  Creat the pursuit flatmap
    final_data = final_mat[pursuit_idx,...]
    param_pursuit = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                     'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                     'cortex_type': 'VertexRGB','description': 'final map',
                     'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                     'cbar_label': '', 'with_labels': True}
    maps_names.append('pursuit')
    
    #  Creat the saccade flatmap
    final_data = final_mat[saccade_idx,...]
    param_saccade = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                     'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                     'cortex_type': 'VertexRGB','description': 'final map',
                     'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                     'cbar_label': '', 'with_labels': True}
    maps_names.append('saccade')
    
    #  Creat the pursuit_and_saccade flatmap
    final_data = final_mat[pursuit_and_saccade_idx,...]
    param_pursuit_and_saccade = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                                  'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                                  'cortex_type': 'VertexRGB','description': 'final map',
                                  'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                                  'cbar_label': '', 'with_labels': True}
    maps_names.append('pursuit_and_saccade')
    
    #  Creat the vision flatmap
    final_data = final_mat[vision_idx,...]
    param_vision = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                    'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                    'cortex_type': 'VertexRGB','description': 'final map',
                    'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                    'cbar_label': '', 'with_labels': True}
    maps_names.append('vision')
    
    #  Creat the vision_and_pursuite flatmap
    final_data = final_mat[vision_and_pursuit_idx,...]
    param_vision_and_pursuit = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                                  'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                                  'cortex_type': 'VertexRGB','description': 'final map',
                                  'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                                  'cbar_label': '', 'with_labels': True}
    maps_names.append('vision_and_pursuit')
    
    #  Creat the vision_and_saccade flatmap
    final_data = final_mat[vision_and_saccade_idx,...]
    param_vision_and_saccade = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                                 'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                                 'cortex_type': 'VertexRGB','description': 'final map',
                                 'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                                 'cbar_label': '', 'with_labels': True}
    maps_names.append('vision_and_saccade')
    
    #  Creat the vision_and_pursuite_and_pursuit flatmap
    final_data = final_mat[vision_and_pursuit_and_saccade_idx,...]
    param_vision_and_pursuit_and_saccade = {'data': final_data, 'cmap': cmap, 'alpha': final_data, 
                                              'vmin': 0, 'vmax': 7, 'cbar': 'stats', 'cmap_steps': cmap_steps,
                                              'cortex_type': 'VertexRGB','description': 'final map',
                                              'curv_brightness': 0.1, 'curv_contrast': 0.25, 'add_roi': save_svg,
                                              'cbar_label': '', 'with_labels': True}
    maps_names.append('vision_and_pursuit_and_saccade')

    
    
    
    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = 'final_stats_{}'.format(maps_name)
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
    
    
