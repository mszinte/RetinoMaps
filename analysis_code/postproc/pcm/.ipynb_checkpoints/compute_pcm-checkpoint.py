"""
-----------------------------------------------------------------------------------------
compute_pcm.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute population cortical magnification and add to derivatives
Note: 
CM is computed using the geodesic distances (mm) of vertices located within a radius on
the flatten surface (see vertex_cm_rad) and restricted by the ROI boundaries
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[4]: model (e.g. css)
-----------------------------------------------------------------------------------------
Output(s):
New brain volume in derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/amblyo_prf/analysis_code/postproc/pcm
2. run python command
>> python compute_pcm.py [main directory] [project name] [subject] [group] [model]
-----------------------------------------------------------------------------------------
Exemple:
python compute_pcm.py /scratch/mszinte/data amblyo_prf sub-01 327 gauss
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import cortex
import importlib
import json
import numpy as np
import os
import sys
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file
import nibabel as nb
import ipdb
deb = ipdb.set_trace

# define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
task = analysis_info['task']
high_pass_type = analysis_info['high_pass_type']
vert_dist_th = analysis_info['vertex_pcm_rad']
rois = analysis_info["rois"]

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
model = sys.argv[5]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# derivatives settings
rsq_idx, ecc_idx, size_idx, x_idx, y_idx = 0, 1, 4, 7, 8

for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    
    # define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, 
                                                           subject, format_)
    fit_dir = "{}/fit".format(prf_dir)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)

    if format_ == 'fsnative':
        deriv_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-deriv_{}_gridfit.func.gii'.format(
            prf_deriv_dir, subject, task, model)
        deriv_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-deriv_{}_gridfit.func.gii'.format(
            prf_deriv_dir, subject, task, model)
        deriv_mat = load_surface_pycortex(L_fn=deriv_avg_fn_L, 
                                          R_fn=deriv_avg_fn_R,)
        
    elif format_ == '170k':
        deriv_avg_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-deriv_{}_gridfit.dtseries.nii'.format(
            prf_deriv_dir, subject, task, model)
        deriv_mat = load_surface_pycortex(brain_fn=deriv_avg_fn)
        save_svg = False

    # parameters
    vert_rsq_data = deriv_mat[rsq_idx, ...]
    vert_x_data = deriv_mat[x_idx, ...]
    vert_y_data = deriv_mat[y_idx, ...]
    vert_size_data = deriv_mat[size_idx, ...]
    vert_ecc_data = deriv_mat[ecc_idx, ...]
    
    # create empty results
    vert_cm = np.zeros(vert_num)*np.nan
    
    # get surfaces for each hemisphere
    surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(subject, "flat")]
    surf_lh, surf_rh = surfs[0], surfs[1]

    # get the vertices number per hemisphere
    lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]
    vert_num = lh_vert_num + rh_vert_num

    # get a dict with the surface vertices contained in each ROI
    roi_verts_dict = cortex.utils.get_roi_verts(subject, mask=False)
    #### TO REPLACE BY YOUR NEW FUNCTION TO GET ROIS FROM NPZ IN CASE OF SUB-170K

    for roi in rois:
        
        # find ROI vertex
        roi_vert_lh_idx = roi_verts_dict[roi][roi_verts_dict[roi]<lh_vert_num]
        roi_vert_rh_idx = roi_verts_dict[roi][roi_verts_dict[roi]>=lh_vert_num]
        roi_surf_lh_idx = roi_vert_lh_idx
        roi_surf_rh_idx = roi_vert_rh_idx-lh_vert_num

        # get mean distance of surounding vertices included in threshold
        vert_lh_rsq, vert_lh_size = vert_rsq_data[:lh_vert_num], vert_size_data[:lh_vert_num]
        vert_lh_x, vert_lh_y = vert_x_data[:lh_vert_num], vert_y_data[:lh_vert_num]
        vert_rh_rsq, vert_rh_size = vert_rsq_data[lh_vert_num:], vert_size_data[lh_vert_num:]
        vert_rh_x, vert_rh_y = vert_x_data[lh_vert_num:], vert_y_data[lh_vert_num:]

        for hemi in ['lh','rh']:
            if hemi == 'lh':
                surf = surf_lh
                roi_vert_idx, roi_surf_idx = roi_vert_lh_idx, roi_surf_lh_idx
                vert_rsq, vert_x, vert_y, vert_size = vert_lh_rsq, vert_lh_x, vert_lh_y, vert_lh_size
            elif hemi == 'rh':
                surf = surf_rh
                roi_vert_idx, roi_surf_idx = roi_vert_rh_idx, roi_surf_rh_idx
                vert_rsq, vert_x, vert_y, vert_size = vert_rh_rsq, vert_rh_x, vert_rh_y, vert_rh_size

            print('ROI -> {} / Hemisphere -> {}'.format(roi, hemi))

            for i, (vert_idx, surf_idx) in enumerate(zip(roi_vert_idx, roi_surf_idx)):

                if vert_rsq[surf_idx] > 0:

                    # get geodesic distances (mm)
                    try :
                        geo_patch = surf.get_geodesic_patch(radius=vert_dist_th, vertex=surf_idx)
                    except Exception as e:
                        print("Vertex #{}: error: {} within {} mm".format(vert_idx, e, vert_dist_th))
                        geo_patch['vertex_mask'] = np.zeros(surf.pts.shape[0]).astype(bool)
                        geo_patch['geodesic_distance'] = []

                    vert_dist_th_idx  = geo_patch['vertex_mask']
                    vert_dist_th_dist = np.ones_like(vert_dist_th_idx)*np.nan
                    vert_dist_th_dist[vert_dist_th_idx] = geo_patch['geodesic_distance']

                    # exclude vextex out of roi
                    vert_dist_th_not_in_roi_idx = [idx for idx in np.where(vert_dist_th_idx)[0] if idx not in roi_surf_idx]
                    vert_dist_th_idx[vert_dist_th_not_in_roi_idx] = False
                    vert_dist_th_dist[vert_dist_th_not_in_roi_idx] = np.nan

                    if np.sum(vert_dist_th_idx) > 0:

                        # compute average geodesic distance excluding distance to itself (see [1:])
                        vert_geo_dist_avg = np.nanmean(vert_dist_th_dist[1:])

                        # get prf parameters of vertices in geodesic distance threshold
                        vert_ctr_x, vert_ctr_y = vert_x[surf_idx], vert_y[surf_idx]
                        vert_dist_th_idx[surf_idx] = False
                        vert_srd_x, vert_srd_y = np.nanmean(vert_x[vert_dist_th_idx]), np.nanmean(vert_y[vert_dist_th_idx])

                        # compute prf center suround distance (deg)
                        vert_prf_dist = np.sqrt((vert_ctr_x - vert_srd_x)**2 + (vert_ctr_y - vert_srd_y)**2)

                        # compute cortical magnification in mm/deg (surface distance / pRF positon distance)
                        vert_cm[vert_idx] = vert_geo_dist_avg/vert_prf_dist
    
    
    # save as last entry in pRF derivative file
    deriv_mat_new = np.zeros((deriv_mat.shape[0]+1, deriv_mat.shape[1]))*np.nan
    deriv_mat_new[0:-1,...] = deriv_mat
    deriv_mat_new[-1,...] = vert_cm
    ### TO SAVE WITH NEW save_surface_pycortex WHICH WILL SPLIT DATA FOR GIFTI BUT NOT CIFTI

# Define permission cmd
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(main_dir, project_dir, group))