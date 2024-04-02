"""
-----------------------------------------------------------------------------------------
add_vertex_surface_to_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
add the surface of each vertex on the final tsv 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# sub-all tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
2. run python command
python finals_figures.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
python add_vertex_surface_to_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import cortex
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import get_rois, calculate_vertex_areas

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
rois = analysis_info['rois']


for format_ in formats : 
    print(format_)
    if format_ == 'fsnative':
        pycortex_subject = subject
        tsv_suffix_list = ['derivatives']
        atlas_name = None
        surf_size = None
        
    elif format_ == '170k':
        pycortex_subject = 'sub-170k'
        tsv_suffix_list = ['derivatives_group', 'derivatives']
        atlas_name = 'mmp_group'
        surf_size = '59k'
        
    # Get surfaces     
    surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(pycortex_subject, "flat")]
    surf_lh, surf_rh = surfs[0], surfs[1]
    lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]
    
    for tsv_suffix in tsv_suffix_list :
        print(tsv_suffix)
        # Get tsv 
        prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv/{}_css-prf_{}.tsv'.format(main_dir, project_dir, subject, format_, subject, tsv_suffix)
        tsv_df = pd.read_table(prf_tsv_fn)
        tsv_df['vertex_surf'] = np.nan


        if format_ == '170k' and tsv_suffix == 'derivatives':
            rois = analysis_info['mmp_rois']
            atlas_name = 'mmp'
        else:
            rois = analysis_info['rois']
            
        # Get rois mask 
        roi_verts_dict_L, roi_verts_dict_R = get_rois(pycortex_subject, return_concat_hemis=False, return_hemi=False, rois=rois, mask=True, atlas_name=atlas_name, surf_size=surf_size)

        for hemi in ['hemi-L', 'hemi-R']:
            print(hemi)
            if hemi == 'hemi-L': 
                surf = surf_lh
                roi_verts_dict = roi_verts_dict_L
                
            elif hemi == 'hemi-R':
                surf = surf_rh
                roi_verts_dict = roi_verts_dict_R
                
            vert_surf = calculate_vertex_areas(surface=surf, mask=None)

            for roi in roi_verts_dict_L.keys():
                print(roi)
                vert_surf_roi = vert_surf[roi_verts_dict[roi]]
                tsv_df.loc[(tsv_df['rois'] == roi) & (tsv_df['hemi'] == hemi), 'vertex_surf'] = vert_surf_roi


                                      
        print(prf_tsv_fn, 'is done')    
        tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
        
# Define permission cmd
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(main_dir, project_dir, group))