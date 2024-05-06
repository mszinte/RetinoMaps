# RetinoMaps
## About
---
We study cerebral areas implicated in both vision and eye movements.</br>
All the analysis are down in surface with both fsnative ad 170k space.
This repository contain all code allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>

---
## Authors (alphabetic order): 
---
Uriel LASCOMBES, Sina KLING, Guillaume MASSON & Martin SZINTE

## Main dependencies
---
[dcm2niix](https://github.com/rordenlab/dcm2niix); 
[PyDeface](https://github.com/poldracklab/pydeface); 
[fMRIprep](https://fmriprep.org/en/stable/); 
[pRFpy](https://github.com/VU-Cog-Sci/prfpy); 
[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/);
[FFmpeg](https://ffmpeg.org/);
[FSL](https://fsl.fmrib.ox.ac.uk);
[Inkscape](https://inkscape.org/);
[workbench](https://humanconnectome.org/software/connectome-workbench)
</br>


## **Data analysis**
---

### To Do 


## Pre-processing
---
#### BIDS
- [x] Copy relevant data from PredictEye [copy_data.py](analysis_code/preproc/bids/bids_copy_data.sh) 
- [x] Change the 'task' to 'task_condition' coulumn name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](analysis_code/preproc/bids/correct_events_files.ipynb)
- [x] Deface participants t1w image [deface_sbatch.py](analysis_code/preproc/bids/deface_sbatch.py) 
    </br>Note: run script for each subject separately.
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]
    </br>Note: for the webpage, use FireFox and wait for at least 30 min, even if nothing seems to happen.

#### Structural preprocessing
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Create sagital view video before manual edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] Manual edit of brain segmentation [pial_edits.sh](analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] Create sagital view video before after edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] Make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [x] Flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, run correlations, average and leave-one-out average [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py) 
- [x] Average inter-run correlations of each subject in 170k template [170k_corr_averaging.py](analysis_code/preproc/functional/170k_corr_averaging.py) 
- [x] Make timeseries inter-run correlations flatmaps with pycortex [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py)

## Post-processing
---
### PRF analysis

#### Gaussian fit
- [x] Create the visual matrix design [vdm_builder.ipynb](analysis_code/postproc/prf/fit/vdm_builder.ipynb)
- [x] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py) 
- [x] Average pRF derivatives from all subjects in 170k template [170k_averaging.py](analysis_code/postproc/prf/postfit/170k_averaging.py) 
- [x] Make pRF derivatives maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py)

#### Rois
- [x] Draw ROIs on individual fsnative using Inkscape
- [x] Import ROIs from the MMP atlas for the hcp 170k template [create_hcp_rois.ipynb](analysis_code/atlas/dev/create_hcp_rois.ipynb)
- [x] Create 170k MMP rois masks [create_mmp_rois_atlas.py](analysis_code/atlas/create_mmp_rois_atlas.py)
- [x] Make ROIS files [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

#### Css fit
- [x] CSS fit within the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [compute_css_stats.py](analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification [css_pcm_sbatch.py](analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make CSS fit derivatives and pcm maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [x] Make TSV with CSS fit derivatives, pcm and statistics [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make pRF derivatives and pcm main figures [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Make subject WEBGL with pycortex [pycortex_webgl_css.py](analysis_code/postproc/prf/webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](disks/meso_H/projects/amblyo_prf/analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](analysis_code/postproc/prf/webgl/publish_webgl.py)

#### For later
- [x] Average pRF derivatives from all subjects in 170k template [170k_averaging.py](analysis_code/postproc/prf/postfit/170k_averaging.py) 

### **Final pycortex maps** 
- [x] Make pRF derivatives and pcm maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py)

### **GLM analysis**
- [x] Run Glm for the differents tasks [glm_sbatch.py](analysis_code/postproc/glm/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](analysis_code/postproc/glm/postfit/compute_glm_stats.py)

### **Statistical analysis**
- [x] Make final statistique image [stats_final.py](analysis_code/postproc/stats/stats_final.py)
- [x] Average statisticales derivatives from all subjects in 170k template [170k_stats_averaging.py](analysis_code/postproc/stats/170k_stats_averaging.py) 
- [x] Make final satistiques maps maps with pycortex [pycortex_maps_stats_final.py](analysis_code/postproc/stats/pycortex_maps_stats_final.py)

### **Group analysis**
- [x] Merge all subjects pycortex maps [pdf_maps.py](analysis_code/postproc/group_analysis/pdf_maps.py)
- [x] Make group tsv [make_group_tsv.py](analysis_code/postproc/group_analysis/make_group_tsv.py)

### **Final ploting**
- [x] Make figures for all subjects [finals_figures.py](analysis_code/postproc/result_analysis/finals_figures.py)


## Utils 



## On working
---
 
### webgl
- [ ] Make webgl with the pycortex dataset [pycortex_maps.py](analysis_code/postproc/prf/webgl/pycortex_webgl.py)
- [ ] Send the files [send_data.sh](analysis_code/postproc/prf/webgl/send_data.sh)

