"""
-----------------------------------------------------------------------------------------
make_rois_img.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make Cifti and Gifti object with rois 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Combined estimate nifti file and pRF derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit
2. run python command
>> python make_rois_img.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
python make_rois_img.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import numpy as np
import nibabel as nb

# personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import make_surface_image, load_surface
from pycortex_utils import get_rois

with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
rois = analysis_info["rois"]
formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
subjects = analysis_info['subjects']

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]


for subject in subjects :
    print(subject)
    for format_, extension in zip(formats, extensions): 
        print(format_)
        rois_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(main_dir, project_dir, subject, format_)
        os.makedirs(rois_dir, exist_ok=True)
        
        if format_ == 'fsnative':
            # Load rois 
            roi_verts_dict_L, roi_verts_dict_R = get_rois(subject, 
                                              return_concat_hemis=False, 
                                              rois=rois, 
                                              mask=True, 
                                              atlas_name=None, 
                                              surf_size=None)
            for hemi in ['hemi-L','hemi-R']:
                if hemi == 'hemi-L':roi_verts_dict = roi_verts_dict_L
                elif hemi == 'hemi-R':roi_verts_dict = roi_verts_dict_R
                array_rois = np.zeros(len(next(iter(roi_verts_dict.values()))), dtype=int)  
                for i, (key, mask) in enumerate(roi_verts_dict.items(), 1):
                    array_rois[mask] = i
                    
                # Load data to have source img
                data_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/fit'.format(main_dir, 
                                                                            project_dir, 
                                                                            subject, 
                                                                            format_)


                
                data_fn = '{}_task-{}_{}_fmriprep_dct_avg_prf-fit_gauss_gridfit.{}'.format(subject, prf_task_name, hemi, extension)
                img, data = load_surface(fn='{}/{}'.format(data_dir, data_fn))
                
                # Exporte data 
                rois_fn = '{}_{}_rois.{}'.format(subject, hemi, extension)
                array_rois = array_rois.reshape(1, -1)
                rois_img = make_surface_image(data=array_rois, source_img=img, maps_names=['rois'])
                nb.save(rois_img, '{}/{}'.format(rois_dir, rois_fn))
   
                
                

                
        elif format_ == '170k':
            roi_verts_dict = get_rois(subject, 
                                      return_concat_hemis=True, 
                                      rois=rois, 
                                      mask=True, 
                                      atlas_name='mmp_group', 
                                      surf_size='170k')

            array_rois = np.zeros(len(next(iter(roi_verts_dict.values()))), dtype=int)  
            for i, (key, mask) in enumerate(roi_verts_dict.items(), 1):
                array_rois[mask] = i

            # Load data to have source img
            data_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/fit'.format(main_dir, 
                                                                        project_dir, 
                                                                        subject, 
                                                                        format_)


            
            data_fn = '{}_task-{}_fmriprep_dct_avg_prf-fit_gauss_gridfit.{}'.format(subject, prf_task_name, extension)
            img, data = load_surface(fn='{}/{}'.format(data_dir, data_fn))
            
            # Exporte data 
            rois_fn = '{}_rois.{}'.format(subject, extension)
            array_rois = array_rois.reshape(1, -1)
            rois_img = make_surface_image(data=array_rois, source_img=img, maps_names=['rois'])
            nb.save(rois_img, '{}/{}'.format(rois_dir, rois_fn))

# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))