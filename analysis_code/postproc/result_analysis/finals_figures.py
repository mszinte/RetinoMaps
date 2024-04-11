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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/result_analysis
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
sys.path.append("{}/../../utils".format(os.getcwd()))
from plot_utils import prf_violins_plot, prf_ecc_size_plot, prf_polar_plot, prf_contralaterality_plot, prf_ecc_pcm_plot, surface_rois_categories_plot ,categories_proportions_roi_plot, surface_rois_all_categories_plot

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Debug
import ipdb
deb = ipdb.set_trace

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
subjects  += ['group', 'sub-170k']
# subjects  += ['sub-170k']


formats = analysis_info['formats']
extensions = analysis_info['extensions']




# Settings 
ecc_th = [0, 20]
size_th= [0.1, 20]
rsq_th = [0, 1]
pcm_th = [0, 20]

for subject in subjects : 
    if subject == 'sub-170k':
        formats = ['170k']
        extensions = ['dtseries.nii']
    for format_, extension in zip(formats, extensions):
        print('making {} figures'.format(subject))
        if format_ == 'fsnative':
            tsv_suffix = 'derivatives'
        elif format_ == '170k':
            tsv_suffix = 'derivatives_group'
            
            
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(main_dir, 
                                                                    project_dir, 
                                                                    subject, 
                                                                    format_)
        
        fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(main_dir, 
                                                                        project_dir, 
                                                                        subject, 
                                                                        format_)
        os.makedirs(fig_dir, exist_ok=True)
        
        # make figures 
        data = pd.read_table('{}/{}_css-prf_{}.tsv'.format(tsv_dir,subject, tsv_suffix))
        fig1 = prf_violins_plot(data, subject,fig_height=1080, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th, pcm_th=pcm_th)
        fig2 = prf_ecc_size_plot(data, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
        figures, hemis = prf_polar_plot(data, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
        fig3 = prf_contralaterality_plot(data, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
        fig4 = prf_ecc_pcm_plot(data, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, pcm_th=pcm_th, rsq_th=rsq_th)
        fig5 = surface_rois_all_categories_plot(data, subject, fig_height=1080, fig_width=1920)
        fig6 = surface_rois_categories_plot(data, subject, fig_height=1080, fig_width=1920)
        fig7 = categories_proportions_roi_plot(data, subject, fig_height=300, fig_width=1920)
        
        #save figures 
        fig1.write_image("{}/{}_prf_rsq_size_n_pcm.pdf".format(fig_dir, subject))
        fig2.write_image("{}/{}_prf_size_ecc.pdf".format(fig_dir, subject)) 
        fig3.write_image("{}/{}_contralaterality.pdf".format(fig_dir, subject)) 
        fig4.write_image("{}/{}_prf_pcm_ecc.pdf".format(fig_dir, subject)) 
        fig5.write_image("{}/{}_surface_rois_all_categories.pdf".format(fig_dir, subject)) 
        fig6.write_image("{}/{}_surface_rois_categories.pdf".format(fig_dir, subject)) 
        fig7.write_image("{}/{}_categories_proportions_roi.pdf".format(fig_dir, subject)) 
        
        for i, (figure, hemi) in enumerate(zip(figures, hemis), start=1):
    
            figure.write_image("{}/{}_subplot_polar_{}.pdf".format(fig_dir, subject, hemi))

categories_to_plot = ['vision',  'vision_and_pursuit_and_saccade','pursuit_and_saccade']

for subject in subjects : 
    if subject == 'sub-170k':
        formats = ['170k']
        extensions = ['dtseries.nii']
    for format_, extension in zip(formats, extensions):
        
        print('making {} figures'.format(subject))
        if format_ == 'fsnative':
            tsv_suffix = 'derivatives'
        elif format_ == '170k':
            tsv_suffix = 'derivatives_group'
            
            
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(main_dir, 
                                                                   project_dir, 
                                                                   subject, 
                                                                   format_)
        
        fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(main_dir, 
                                                                       project_dir, 
                                                                       subject, 
                                                                       format_)
        os.makedirs(fig_dir, exist_ok=True)
        
        # make figures 
        data = pd.read_table('{}/{}_css-prf_{}.tsv'.format(tsv_dir,subject, tsv_suffix))
        for categorie_to_plot in categories_to_plot : 
            print(categorie_to_plot)
            fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures/{}'.format(main_dir, 
                                                                              project_dir, 
                                                                              subject, 
                                                                              format_, 
                                                                              categorie_to_plot)
            os.makedirs(fig_dir, exist_ok=True)
            df_categorie = data.loc[data.stats_final == categorie_to_plot]
            
            # fig1 = prf_violins_plot(df_categorie, subject,fig_height=1080, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th, pcm_th=pcm_th)

            # fig2 = prf_ecc_size_plot(df_categorie, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
            # figures, hemis = prf_polar_plot(df_categorie, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
            # fig3 = prf_contralaterality_plot(df_categorie, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
            # fig4 = prf_ecc_pcm_plot(df_categorie, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, pcm_th=pcm_th, rsq_th=rsq_th)

            
            # #save figures 
            # fig1.write_image("{}/{}_prf_rsq_size_n_pcm_{}.pdf".format(fig_dir, subject, categorie_to_plot))
            # fig2.write_image("{}/{}_prf_size_ecc_{}.pdf".format(fig_dir, subject, categorie_to_plot)) 
            # fig3.write_image("{}/{}_contralaterality_{}.pdf".format(fig_dir, subject, categorie_to_plot)) 
            # fig4.write_image("{}/{}_prf_pcm_ecc_{}.pdf".format(fig_dir, subject, categorie_to_plot)) 

            
            # for i, (figure, hemi) in enumerate(zip(figures, hemis), start=1):
        
            #     figure.write_image("{}/{}_subplot_polar_{}_{}.pdf".format(fig_dir, subject, hemi, categorie_to_plot))
    
            try:
                fig1 = prf_violins_plot(df_categorie, subject,fig_height=1080, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th, pcm_th=pcm_th)
                fig1.write_image("{}/{}_prf_rsq_size_n_pcm_{}.pdf".format(fig_dir, subject, categorie_to_plot))
            except Exception:
                pass
            
            try:
                fig2 = prf_ecc_size_plot(df_categorie, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
                fig2.write_image("{}/{}_prf_size_ecc_{}.pdf".format(fig_dir, subject, categorie_to_plot))
            except Exception:
                pass
            
            try:
                fig3 = prf_contralaterality_plot(df_categorie, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
                fig3.write_image("{}/{}_contralaterality_{}.pdf".format(fig_dir, subject, categorie_to_plot))
            except Exception:
                pass
            
            try:
                fig4 = prf_ecc_pcm_plot(df_categorie, subject, fig_height=400, fig_width=800, ecc_th=ecc_th, pcm_th=pcm_th, rsq_th=rsq_th)
                fig4.write_image("{}/{}_prf_pcm_ecc_{}.pdf".format(fig_dir, subject, categorie_to_plot))
            except Exception:
                pass
            
            try:
                figures, hemis = prf_polar_plot(df_categorie, subject, fig_height=300, fig_width=1920, ecc_th=ecc_th, size_th=size_th, rsq_th=rsq_th)
                for i, (figure, hemi) in enumerate(zip(figures, hemis), start=1):
                    figure.write_image("{}/{}_subplot_polar_{}_{}.pdf".format(fig_dir, subject, hemi, categorie_to_plot))
            except Exception:
                pass


        
# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))
    
    
    
    
    
    
    
    
    
    