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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/pcm
2. run python command
>> python compute_pcm.py [main directory] [project name] [subject] [group] [model]
-----------------------------------------------------------------------------------------
Exemple:
python compute_pcm.py /scratch/mszinte/data RetinoMaps sub-01 327 css
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# general imports
import os
import sys
import json
import glob
import cortex
import importlib
import numpy as np
import pandas as pd
import nibabel as nb

# Personal iports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file, load_surface_pycortex, get_rois, make_image_pycortex

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
model = sys.argv[5]

if model == 'gauss':
    model = 'gauss_gridfit'
elif model == 'css':
    model = 'avg_css'

# define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info["task_names"]
task = tasks[2]

vert_dist_th = analysis_info['vertex_pcm_rad']
formats = analysis_info['formats']

# # Set pycortex db and colormaps
# cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
# set_pycortex_config_file(cortex_dir)
# importlib.reload(cortex)

# derivatives settings
rsq_idx, ecc_idx, size_idx, x_idx, y_idx = 0, 1, 4, 7, 8

for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):

    # define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject, format_)
    fit_dir = "{}/fit".format(prf_dir)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    

    if format_ == 'fsnative':
        rois = analysis_info['rois']
        atlas_name = None 
        surf_size = None        
        deriv_avg_fn_L = glob.glob('{}/{}*hemi-L*prf-deriv-*{}.func.gii'.format(prf_deriv_dir, 
                                                                                   subject, 
                                                                                   model))
        
        deriv_avg_fn_R = glob.glob('{}/{}*hemi-R*prf-deriv-*{}.func.gii'.format(prf_deriv_dir, 
                                                                                   subject, 
                                                                                   model))
   
        
        results = load_surface_pycortex(L_fn=deriv_avg_fn_L[0], 
                                        R_fn=deriv_avg_fn_R[0], 
                                        return_img=True)
        
        deriv_mat = results['data_concat'] 
        img_L = results['img_L'] 
        img_R = results['img_R']
        
    elif format_ == '170k':
        rois = analysis_info['mmp_rois']
        atlas_name = 'mmp'
        surf_size = '59k'
        deriv_avg_fn = glob.glob('{}/{}*prf-deriv-*{}.dtseries.nii'.format(prf_deriv_dir, 
                                                                     subject, 
                                                                     model)) 
        results = load_surface_pycortex(brain_fn=deriv_avg_fn[0],
                                        return_img=True,
                                        return_59k_mask=True,  
                                        return_source_data=True)
        
        deriv_mat = results['data_concat']
        mask_59k = results['mask_59k']
        deriv_mat_170k = results['source_data'] 
        img = results['img']

    # get surfaces for each hemisphere
    surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(pycortex_subject, 
                                                                      "flat")]
    surf_lh, surf_rh = surfs[0], surfs[1]
    # get the vertices number per hemisphere
    lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]
    vert_num = lh_vert_num + rh_vert_num
    
    # get a dicst with the surface vertices contained in each ROI
    roi_verts_dict = get_rois(pycortex_subject, 
                              return_concat_hemis=True, 
                              rois=rois, 
                              mask=False, 
                              atlas_name=atlas_name, 
                              surf_size=surf_size)

    
    # derivatives settings
    rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
        amp_idx, baseline_idx, x_idx, y_idx, hrf_1_idx, hrf_2_idx = \
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    if model == 'gauss':
        loo_rsq_idx = 11
    elif model == 'dn':
        srf_amp_idx = 11
        sff_size = 12
        neural_baseline_idx = 13
        surround_baseline = 14
        loo_rsq_idx = 15
    elif model == 'css':
        n_idx = 11
        loo_rsq_idx = 12
    
    # parameters
    vert_rsq_data = deriv_mat[rsq_idx, ...]
    vert_x_data = deriv_mat[x_idx, ...]
    vert_y_data = deriv_mat[y_idx, ...]
    vert_size_data = deriv_mat[size_idx, ...]
    vert_ecc_data = deriv_mat[ecc_idx, ...]
    
    # create empty results
    vert_cm = np.zeros(vert_num)*np.nan
    
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
    
            desc = 'ROI -> {} / Hemisphere -> {}'.format(roi, hemi)
            print(desc)
            for i, (vert_idx, surf_idx) in enumerate(zip(roi_vert_idx, roi_surf_idx)):

                if vert_rsq[surf_idx] > 0:
    
                    # get geodesic distances (mm)
                    try :
                        geo_patch = surf.get_geodesic_patch(radius=vert_dist_th, 
                                                            vertex=surf_idx)
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

    deriv_mat_new = np.zeros((deriv_mat.shape[0]+1, deriv_mat.shape[1]))*np.nan
    deriv_mat_new[0:-1,...] = deriv_mat
    deriv_mat_new[-1,...] = vert_cm
    
    # Exporte Data
    if model == 'gauss_gridfit':
            maps_names = ['rsq', 'ecc', 'polar_real', 'polar_imag', 'size', 
                          'amplitude', 'baseline', 'x', 'y', 'hrf_1',' hrf_2', 
                          'pcm']
    elif model == 'avg_css':
        maps_names = ['prf_rsq', 'prf_ecc', 'polar_real', 'polar_imag', 
                      'prf_size', 'amplitude', 'baseline', 'prf_x','prf_y',
                      ' hrf_1', 'hrf_2','prf_n', 'prf_loo_r2', 'pcm']

    if format_ == 'fsnative':

        new_img_L, new_img_R  = make_image_pycortex(data=deriv_mat_new, 
                                                    maps_names=maps_names,
                                                    img_L=img_L, 
                                                    img_R=img_R, 
                                                    lh_vert_num=lh_vert_num, 
                                                    rh_vert_num=rh_vert_num, 
                                                    img=None, 
                                                    brain_mask_59k=None)
        
        deriv_avg_fn_L = deriv_avg_fn_L[0].split('/')[-1].replace('deriv', 'derivs_pcm')  
        deriv_avg_fn_R = deriv_avg_fn_R[0].split('/')[-1].replace('deriv', 'derivs_pcm')
        nb.save(new_img_L, '{}/{}'.format(prf_deriv_dir,deriv_avg_fn_L))
        nb.save(new_img_R, '{}/{}'.format(prf_deriv_dir,deriv_avg_fn_R))
        
        
        prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv/{}_css-prf_derivatives.tsv'.format(main_dir, 
                                                                                                 project_dir, 
                                                                                                 subject, 
                                                                                                 format_, 
                                                                                                 subject)
        tsv_df = pd.read_table(prf_tsv_fn)
        
        
        tsv_df['pcm'] = np.nan
    
        for roi in roi_verts_dict.keys():
            data_roi = deriv_mat_new[-1, roi_verts_dict[roi]]
            tsv_df.loc[tsv_df['rois'] == roi, 'pcm'] = data_roi
        
        tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
                
    elif format_ == '170k':
        new_img = make_image_pycortex(data=deriv_mat_new, 
                              maps_names=None,
                              img_L=None, 
                              img_R=None, 
                              lh_vert_num=None, 
                              rh_vert_num=None, 
                              img=img, 
                              brain_mask_59k=mask_59k)
        deriv_avg_fn = deriv_avg_fn[0].split('/')[-1].replace('deriv', 'derivs_pcm')

        nb.save(new_img, '{}/{}'.format(prf_deriv_dir,deriv_avg_fn))
        
        concat_rois_list = [analysis_info['mmp_rois'], analysis_info['rois']]
        for n_list, rois_list in enumerate(concat_rois_list):
            rois = rois_list
            if 'LO' in rois_list:
                atlas_name = 'mmp_group'
                tsv_suffix = 'derivatives_group'
            else:
                atlas_name = 'mmp'
                tsv_suffix = 'derivatives'
                
            roi_verts_dict = get_rois(pycortex_subject, 
                                      return_concat_hemis=True, 
                                      rois=rois, 
                                      mask=False, 
                                      atlas_name=atlas_name, 
                                      surf_size=surf_size)
        
            prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv/{}_css-prf_{}.tsv'.format(main_dir, 
                                                                                            project_dir, 
                                                                                            subject, 
                                                                                            format_, 
                                                                                            subject, 
                                                                                            tsv_suffix)
            tsv_df = pd.read_table(prf_tsv_fn)
            tsv_df['pcm'] = np.nan
    
            for roi in roi_verts_dict.keys():
                data_roi = deriv_mat_new[-1, roi_verts_dict[roi]]
                tsv_df.loc[tsv_df['rois'] == roi, 'pcm'] = data_roi
            
            tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(main_dir, project_dir, group))