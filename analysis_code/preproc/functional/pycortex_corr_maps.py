"""
-----------------------------------------------------------------------------------------
pycortex_maps.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save map in svg (y/n)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/preproc/functional/
2. run python command
>> python pycortex_corr_maps.py [main directory] [project name] [subject num] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
python pycortex_corr_maps.py ~/disks/meso_shared RetinoMaps sub-02 n
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
import nibabel as nb
import numpy as np
import os
import sys
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file
deb = ipdb.set_trace

# Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
xfm_name = analysis_info["xfm_name"]

high_pass_type = analysis_info['high_pass_type']
# tasks = analysis_info['task_names']
tasks = ['PurVELoc','SacVELoc']


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
 
# Define directories and fn
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
corr_dir = "{}/{}/derivatives/pp_data/{}/func/fmriprep_dct_corr/fsnative".format(main_dir, project_dir, subject)
flatmaps_corr_dir = '{}/{}/derivatives/pp_data/{}/corr/pycortex/flatmaps_corr'.format(main_dir, project_dir, subject)
datasets_corr_dir = '{}/{}/derivatives/pp_data/{}/corr/pycortex/datasets_corr'.format(main_dir, project_dir, subject)




os.makedirs(flatmaps_corr_dir, exist_ok=True)
os.makedirs(datasets_corr_dir, exist_ok=True)

flatmaps_dir = flatmaps_corr_dir
datasets_dir = datasets_corr_dir

deriv_fn_label = 'corr'

# Set pycortex db and colormaps
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# Maps settings



 
cmap_corr = 'RdBu_r_alpha'


# plot scales
corr_scale = [-1, 1]

for task in tasks : 
    print(task)
    deriv_corr_fn_L = "{}/{}_task-{}_hemi-L_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,high_pass_type)
    deriv_corr_fn_R = "{}/{}_task-{}_hemi-R_fmriprep_{}_correlations_bold.func.gii".format(corr_dir, subject, task,high_pass_type)
    

    
    
    print('Creating flatmaps...')
    
            
    maps_names = []
    
    # load data
    img_corr_L = nb.load(deriv_corr_fn_L)
    corr_data_L = [x.data for x in img_corr_L.darrays]
    corr_data_L = np.vstack(corr_data_L) 
    
    img_corr_R = nb.load(deriv_corr_fn_R)
    corr_data_R = [x.data for x in img_corr_R.darrays]
    corr_data_R = np.vstack(corr_data_R) 
    
    corr_data_concat = np.concatenate((corr_data_L,corr_data_R))
    corr_data_vect = corr_data_concat.reshape(corr_data_concat.shape[0])
    
    print('data are load')
    # correlation

    alpha_range = analysis_info["alpha_range"]
    alpha = (corr_data_vect - alpha_range[0])/(alpha_range[1]-alpha_range[0])
    alpha[alpha>1]=1
    alpha = alpha.astype(np.uint8)
    
    
    
    param_correlations = {'data': corr_data_vect, 'cmap': cmap_corr ,'alpha' : alpha,
                  'vmin': corr_scale[0], 'vmax': corr_scale[1], 'cbar': 'discrete', 
                  'cortex_type': 'VertexRGB', 'description': '{} {} correlation'.format(subject,task), 
                  'curv_brightness': 1, 'curv_contrast': 0.1, 'add_roi': save_svg, 'cbar_label': 'Pearson coefficient',
                  'with_labels': True}
    maps_names.append('correlations')
    
    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = '{}_{}'.format(task, maps_name)
        roi_param = {'subject': subject, 'xfmname': None, 'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name = maps_name))
        exec("plt.savefig('{}/{}_task-{}_{}_{}.pdf')".format(flatmaps_dir, subject, task,  maps_name, deriv_fn_label))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_task-{}_{}.hdf".format(datasets_dir, subject, task, deriv_fn_label)
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)
    
    
    
