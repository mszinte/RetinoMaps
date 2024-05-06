"""
-----------------------------------------------------------------------------------------
make_rois_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs-based CSS figures
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import json
import pandas as pd
import ipdb
deb = ipdb.set_trace

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']

# Threshold settings
ecc_th = analysis_info['ecc_th']
size_th = analysis_info['size_th']
rsqr_th = analysis_info['rsqr_th']
pcm_th = analysis_info['pcm_th']
amplitude_th = analysis_info['amplitude_th']
stats_th = analysis_info['stats_th']
if stats_th == 0.05: stats_col = 'corr_pvalue_5pt'
elif stats_th == 0.01: stats_col = 'corr_pvalue_1pt'

# Figure settings
colormap_dict = {'V1': (243, 231, 155),
                 'V2': (250, 196, 132),
                 'V3': (248, 160, 126),
                 'V3AB': (235, 127, 134),
                 'LO': (150, 0, 90), 
                 'VO': (0, 0, 200),
                 'hMT+': (0, 25, 255),
                 'iIPS': (0, 152, 255),
                 'sIPS': (44, 255, 150),
                 'iPCS': (151, 255, 0),
                 'sPCS': (255, 234, 0),
                 'mPCS': (255, 111, 0)
                }
roi_colors = ['rgb({},{},{})'.format(*rgb) for rgb in colormap_dict.values()]
plot_groups = [['V1', 'V2', 'V3'], ['V3AB', 'LO', 'VO'], ['hMT+', 'iIPS', 'sIPS'], ['iPCS', 'sPCS', 'mPCS']]
num_ecc_size_bins = 8
num_ecc_pcm_bins = 8
num_polar_bins = 12
max_ecc = 15
fig_width = 1200

# Format loop
for format_, extension in zip(formats, extensions):
    
    # Create folders and fns
    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
        main_dir, project_dir, subject, format_)
    fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
        main_dir, project_dir, subject, format_)
    os.makedirs(fig_dir, exist_ok=True)

    # Load data
    tsv_fn = '{}/{}_css-all_derivatives.tsv'.format(tsv_dir, subject)
    data = pd.read_table(tsv_fn, sep="\t")
    
    # keep a raw data df 
    data_raw = data.copy()

    # Threshold data (replace by nan)
    data.loc[(data.amplitude < amplitude_th) |                                    # amplitude 
             (data.prf_ecc < ecc_th[0]) | (data.prf_ecc > ecc_th[1]) |            # eccentricity 
             (data.prf_size < size_th[0]) | (data.prf_size > size_th[1]) |        # size
             (data.prf_loo_r2 < rsqr_th) |                                        # loo rsqr
             (data[stats_col] > stats_th)                                         # stats
              ] = np.nan
    data = data.dropna()

    # Stats plot
    fig_fn = "{}/{}_prf_roi_area.pdf".format(fig_dir, subject)
    print('Saving {}'.format(fig_fn))
    fig = prf_roi_area(data=data_raw, fig_width=fig_width, fig_height=600, roi_colors=roi_colors)
    fig.write_image(fig_fn)
    
    # Violins plots
    fig_fn = "{}/{}_prf_violins.pdf".format(fig_dir, subject)
    print('Saving {}'.format(fig_fn))
    fig = prf_violins_plot(data=data, fig_width=fig_width, fig_height=600, rois=rois, roi_colors=roi_colors)
    fig.write_image(fig_fn)

    # Ecc.size plots
    fig_fn = "{}/{}_prf_ecc_size.pdf".format(fig_dir, subject)
    print('Saving {}'.format(fig_fn))
    fig = prf_ecc_size_plot(data=data, fig_width=fig_width, fig_height=400, rois=rois, roi_colors=roi_colors,
                            plot_groups=plot_groups, num_bins=num_ecc_size_bins, max_ecc=max_ecc)
    fig.write_image(fig_fn)

    # Ecc.pCM plot
    data_pcm = data.copy()
    data_pcm.loc[(data_pcm.pcm < pcm_th[0]) | (data_pcm.pcm > pcm_th[1])] = np.nan
    data_pcm = data_pcm.dropna()
    fig_fn = "{}/{}_prf_ecc_pcm.pdf".format(fig_dir, subject)
    fig = prf_ecc_pcm_plot(data_pcm, fig_width=fig_width, fig_height=400, rois=rois, roi_colors=roi_colors,
                           plot_groups=plot_groups, num_bins=num_ecc_pcm_bins, max_ecc=max_ecc)
    print('Saving {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Polar angle distributions
    figs, hemis = prf_polar_plot(data, fig_width=fig_width, fig_height=300, rois=rois, roi_colors=roi_colors, num_bins=num_polar_bins)
    for (fig, hemi) in zip(figs, hemis):
        fig_fn = "{}/{}_prf_polar_angle_{}.pdf".format(fig_dir, subject, hemi)
        print('Saving {}'.format(fig_fn))
        fig.write_image(fig_fn)

    # Contralaterality plots
    fig_fn = "{}/{}_contralaterality.pdf".format(fig_dir, subject)
    fig = prf_contralaterality_plot(data, fig_width=fig_width, fig_height=300, rois=rois, roi_colors=roi_colors)
    print('Saving {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Spatial distibution plot
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))