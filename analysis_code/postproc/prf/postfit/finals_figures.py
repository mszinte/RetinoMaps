"""
-----------------------------------------------------------------------------------------
finals_figures.py
-----------------------------------------------------------------------------------------
Goal of the script:
make finals figures for all subjects
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
python finals_figures.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python finals_figures.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports


import os
import sys
import json
import pandas as pd


# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import prf_violins_plot, prf_ecc_size_plot, prf_polar_plot, prf_contralaterality_plot 

# # figure imports
# import plotly.graph_objects as go
# import plotly.express as px

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# general imports
import ipdb
deb = ipdb.set_trace



# load settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
subjects  += ['sub-all']



# Settings 
ecc_th = [0, 15]
size_th= [0.1, 20]
rsq_th = [0.05, 1]

for subject in subjects : 
    print('making {} figures'.format(subject))
    tsv_dir = '{}/{}/derivatives/pp_data/{}/fsnative/prf/tsv'.format(main_dir, 
                                                                     project_dir, 
                                                                     subject)
    
    fig_dir = '{}/{}/derivatives/pp_data/{}/fsnative/prf/figures'.format(main_dir, 
                                                                         project_dir, 
                                                                         subject)
    os.makedirs(fig_dir, exist_ok=True)
    
    data = pd.read_table('{}/{}_task-prf_loo.tsv'.format(tsv_dir,subject))
    fig1 = prf_violins_plot(data, subject, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
    fig2 = prf_ecc_size_plot(data, subject, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
    figures, hemis = prf_polar_plot(data, subject, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
    fig3 = prf_contralaterality_plot(data, subject, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
    
    
    fig1.write_image("{}/{}_prf_rsq_size_n.pdf".format(fig_dir, subject))
    fig2.write_image("{}/{}_prf_size_ecc.pdf".format(fig_dir, subject)) 
    fig3.write_image("{}/{}_contralaterality.pdf".format(fig_dir, subject)) 
    
    for i, (figure, hemi) in enumerate(zip(figures, hemis), start=1):

        figure.write_image("{}/{}_subplot_polar_{}.pdf".format(fig_dir, subject, hemi))
    
# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))
    
    
    
    
    
    
    
    
    
    