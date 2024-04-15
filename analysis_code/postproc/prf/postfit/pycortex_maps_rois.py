"""
-----------------------------------------------------------------------------------------
pycortex_maps_rois.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save in svg (e.g. no)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_rois.py [main directory] [project name] [subject num] [save_in_svg]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_maps_rois.py ~/disks/meso_shared RetinoMaps sub-01 n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug import 
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex

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
if subject == 'sub-170k': save_svg = save_svg
else: save_svg = False

# Define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    rois_dir = "{}/{}/derivatives/pp_data/{}/{}/rois".format(main_dir, project_dir, subject,format_)
    flatmaps_dir = '{}/pycortex/flatmaps_rois'.format(rois_dir)
    datasets_dir = '{}/pycortex/datasets_rois'.format(rois_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    if format_ == 'fsnative':
        roi_fn_L = '{}/{}_hemi-L_rois.func.gii'.format(rois_dir, subject)
        roi_fn_R = '{}/{}_hemi-R_rois.func.gii'.format(rois_dir, subject)
        results = load_surface_pycortex(L_fn=roi_fn_L, 
                                        R_fn=roi_fn_R)
        roi_mat = results['data_concat']
        
    elif format_ == '170k':
        roi_fn = '{}/{}_rois.dtseries.nii'.format(rois_dir, subject)
        results = load_surface_pycortex(brain_fn=roi_fn)
        roi_mat = results['data_concat']
    
    print('Creating flatmaps...')
    maps_names = []

    # rois
    param_rois = {'subject': subject,
                  'data': roi_mat, 
                  'cmap': 'rois_colors', 
                  'alpha': roi_mat, 
                  'vmin': 0, 
                  'vmax': 12, 
                  'cbar': 'rois', 
                  'cortex_type': 'VertexRGB',
                  'description': '',
                  'curv_brightness': 1, 
                  'curv_contrast': 0.1, 
                  'add_roi': save_svg,
                  'cbar_label': '',
                  'with_labels': True}
    maps_names.append('rois')

    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = '{}_{}'.format(prf_task_name, maps_name)
        roi_param = {'subject': pycortex_subject, 
                     'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
        exec("plt.savefig('{}/{}_task-{}_rois.pdf')".format(flatmaps_dir, subject, prf_task_name))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})

    # save dataset
    dataset_file = "{}/{}_task-{}_{}.hdf".format(datasets_dir, subject, prf_task_name, deriv_fn_label)
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)