"""
-----------------------------------------------------------------------------------------
make_tsv_css.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create TSV file with all css analysis output
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
-----------------------------------------------------------------------------------------
Output(s):
TSV file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_stats.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
python compute_css_stats.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
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
import matplotlib.pyplot as plt

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap

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
formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Maps settings 
rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx, \
    size_idx, amp_idx, baseline_idx, x_idx, y_idx, \
    n_idx, loo_rsq_idx, pcm_idx, slope_idx, intercept_idx, \
    rvalue_idx, pvalue_idx, stderr_idx, pvalue_corrected_5pt_idx,  \
    pvalue_corrected_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20

cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

# plot scales
rsq_scale = [0, 1]
ecc_scale = [0, 8]
size_scale = [0, 5]
n_scale = [0, 2]
pcm_scale = [0, 10]

for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    
    # Define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject, format_)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    flatmaps_dir = '{}/pycortex/flatmaps_loo-avg_css'.format(prf_dir)
    datasets_dir = '{}/pycortex/datasets_loo-avg_css'.format(prf_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)

    # Load all data
    if format_ == 'fsnative':
        # Derivatives
        deriv_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_prf-deriv-loo-avg_css.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_prf-deriv-loo-avg_css.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_results = load_surface_pycortex(L_fn=deriv_avg_fn_L, R_fn=deriv_avg_fn_R)
        deriv_mat = deriv_results['data_concat']
        
        # pcm 
        pcm_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_prf-pcm-loo-avg_css.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_prf-pcm-loo-avg_css.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_results = load_surface_pycortex(L_fn=pcm_avg_fn_L, R_fn=pcm_avg_fn_R)
        pcm_mat = pcm_results['data_concat']
        
        # Stats
        stats_avg_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_avg_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_results = load_surface_pycortex(L_fn=stats_avg_fn_L, R_fn=stats_avg_fn_R)
        stats_mat = stats_results['data_concat']
        
    elif format_ == '170k':
        # Derivatives
        deriv_avg_fn = '{}/{}_task-{}_fmriprep_dct_prf-deriv-loo-avg_css.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_results = load_surface_pycortex(brain_fn=deriv_avg_fn)
        deriv_mat = deriv_results['data_concat']

        # pcm
        pcm_avg_fn = '{}/{}_task-{}_fmriprep_dct_prf-pcm-loo-avg_css.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_results = load_surface_pycortex(brain_fn=pcm_avg_fn)
        pcm_mat = pcm_results['data_concat']

        # Stats
        stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_loo-avg_prf-stats.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_results = load_surface_pycortex(brain_fn=stats_avg_fn)
        stats_mat = stats_results['data_concat']
        
    # Combine mat
    all_deriv_mat = np.concatenate((deriv_mat, pcm_mat, stats_mat))
    
    # Threshold mat
    all_deriv_mat_th = all_deriv_mat
    amp_down = all_deriv_mat_th[amp_idx,...] > 0
    rsq_down = all_deriv_mat_th[loo_rsq_idx,...] >= 0
    size_th_down = all_deriv_mat_th[size_idx,...] >= analysis_info['size_th'][0]
    size_th_up = all_deriv_mat_th[size_idx,...] <= analysis_info['size_th'][1]
    ecc_th_down = all_deriv_mat_th[ecc_idx,...] >= analysis_info['ecc_th'][0]
    ecc_th_up = all_deriv_mat_th[ecc_idx,...] <= analysis_info['ecc_th'][1]
    if analysis_info['stats_th'] == 0.05: stats_th_down = all_deriv_mat_th[pvalue_corrected_5pt_idx,...] <= 0.05
    elif analysis_info['stats_th'] == 0.01: stats_th_down = all_deriv_mat_th[pvalue_corrected_1pt_idx,...] <= 0.01
    all_th = np.array((amp_down, 
                       rsq_down,
                       size_th_down,size_th_up, 
                       ecc_th_down, ecc_th_up,
                       stats_th_down)) 
    all_deriv_mat[loo_rsq_idx, np.logical_and.reduce(all_th)==False]=0 # put this to zero to not plot it
# rois_group = analysis_info['rois_group']

# # Make TSV
# for format_, extension in zip(formats, extensions):
#     prf_tsv_dir = "{}/{}/{}/prf/tsv".format(pp_dir, subject, format_)
#     os.makedirs(prf_tsv_dir, exist_ok=True)
#     tsv_fn = '{}/{}_css-prf_derivatives.tsv'.format(prf_tsv_dir, subject)
#     df_rois_brain = pd.DataFrame()
    
#     if format_ == 'fsnative':
#         atlas_name, surf_size = None, None 

#         for hemi in ['hemi-L', 'hemi-R']:
#             brain_data_avg = hemi_data_avg[hemi]
#             roi_verts = get_rois(subject=subject, 
#                                  return_concat_hemis=False, 
#                                  return_hemi=hemi, 
#                                  rois=rois, 
#                                  mask=True, 
#                                  atlas_name=atlas_name, 
#                                  surf_size=surf_size)

#             for roi in roi_verts.keys():
#                 data_dict = {col: brain_data_avg[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
#                 data_dict['roi'] = [roi] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                 data_dict['subject'] = [subject] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                 data_dict['hemi'] = [hemi] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                 df_rois_hemi = pd.DataFrame(data_dict)
#                 df_rois_brain = pd.concat([df_rois_brain, df_rois_hemi], ignore_index=True)

#         df_rois_brain.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)
#         print('Saving tsv: {}'.format(tsv_fn))


        # # Save as tsv
        # prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv/{}_css-prf_derivatives.tsv'.format(
        #     main_dir, project_dir, subject, format_, subject)
        # tsv_df = pd.read_table(prf_tsv_fn)
        # tsv_df['pcm'] = np.nan
    
        # for roi in roi_verts_dict.keys():
        #     data_roi = deriv_mat_new[-1, roi_verts_dict[roi]]
        #     tsv_df.loc[tsv_df['rois'] == roi, 'pcm'] = data_roi
        
        # tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)

#     elif format_ == '170k':
#         atlas_name, surf_size = 'mmp', '170k'
        
#         for rois_group, roi in zip(rois_group, rois):
#             brain_data_avg = hemi_data_avg['170k']
#             roi_verts_L, roi_verts_R = get_rois(subject=subject,
#                                                 return_concat_hemis=False,
#                                                 return_hemi=None,
#                                                 rois=rois_group,
#                                                 mask=True,
#                                                 atlas_name=atlas_name,
#                                                 surf_size=surf_size)
            
#             # combine dict with multiple areas (thanks ChatGPT)
#             roi_verts_L = {roi: np.array([all(value) for value in zip(*roi_verts_L.values())])}
#             roi_verts_R = {roi: np.array([all(value) for value in zip(*roi_verts_R.values())])}
            
#             for hemi in ['hemi-L', 'hemi-R']:
#                 if hemi == 'hemi-L': roi_verts = roi_verts_L
#                 elif hemi == 'hemi-R': roi_verts = roi_verts_R

#                 for roi in roi_verts.keys():
#                     data_dict = {col: brain_data_avg[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
#                     data_dict['roi'] = [roi] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                     data_dict['subject'] = [subject] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                     data_dict['hemi'] = [hemi] * brain_data_avg[:, roi_verts[roi]].shape[1]
#                     df_rois_hemi = pd.DataFrame(data_dict)
#                     df_rois_brain = pd.concat([df_rois_brain, df_rois_hemi], ignore_index=True)
                    
#         df_rois_brain.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)
#         print('Saving tsv: {}'.format(tsv_fn))


            
        # # concat_rois_list = [analysis_info['mmp_rois'], analysis_info['rois']]
        # # for n_list, rois_list in enumerate(concat_rois_list):
        #     # rois = rois_list
        #     # if 'LO' in rois_list:
        #     #     atlas_name = 'mmp_group'
        #     #     tsv_suffix = 'derivatives_group'
        #     # else:
        #     #     atlas_name = 'mmp'
        #     #     tsv_suffix = 'derivatives'
                
        #     roi_verts_dict = get_rois(pycortex_subject, 
        #                               return_concat_hemis=True, 
        #                               rois=rois, 
        #                               mask=False, 
        #                               atlas_name=atlas_name, 
        #                               surf_size=surf_size)
        
        #     prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv/{}_css-prf_{}.tsv'.format(
        #         main_dir, project_dir, subject, format_, subject, tsv_suffix)
            
        #     tsv_df = pd.read_table(prf_tsv_fn)
        #     tsv_df['pcm'] = np.nan
    
        #     for roi in roi_verts_dict.keys():
        #         data_roi = deriv_mat_new[-1, roi_verts_dict[roi]]
        #         tsv_df.loc[tsv_df['rois'] == roi, 'pcm'] = data_roi
            
        #     tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
