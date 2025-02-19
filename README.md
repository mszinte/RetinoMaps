# RetinoMaps
## About
---
*We here study cortical areas involved in both vision and eye movements in 20 healthy controls.</br>*
*All the analyses are done in surface with both **fsnative** ad **HCP 170k** format.</br>*
*This repository contains all codes allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>*

---
## Authors (alphabetic order): 
---
Sina KLING, Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

## To Do 
---
- [ ] take care of glm analysis

## Data analysis
---

### BIDS
- [x] Copy relevant data from PredictEye [copy_data.py](analysis_code/preproc/bids/bids_copy_data.sh) 
- [x] Change the 'task' to 'task_condition' coulumn name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](analysis_code/preproc/bids/correct_events_files.ipynb)
- [x] Deface participants t1w image [deface_sbatch.py](analysis_code/preproc/bids/deface_sbatch.py) 
    </br>Note: run script for each subject separately.
- [x] Validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]
    </br>Note: for the webpage, use FireFox and wait for at least 30 min, even if nothing seems to happen.

### Individual analysis
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

#### Structural preprocessing
- [x] fMRIprep with anat-only option [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] create sagital view video before manual edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] manual edit of brain segmentation [pial_edits.sh](analysis_code/preproc/anatomical/pial_edits.sh)
- [x] FreeSurfer with new brainmask manually edited [freesurfer_pial.py](analysis_code/preproc/anatomical/freesurfer_pial.py)
- [x] create sagital view video before after edit [sagital_view.py](analysis_code/preproc/anatomical/sagital_view.py)
- [x] make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [x] flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [x] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [x] supress bad run [bad_run.py](analysis_code/preproc/functional/bad_run.py)
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, average and leave-one-out average and correlations [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py)
- [x] Compute vertex areas [compute_vertex_area.py](analysis_code/preproc/anatomical/compute_vertex_area.py)

#### Functional postprocessing
Analyses are run on individual participant (**sub-0X**) surface (**fsnative**) or their projection on the HCP cifti format (**170k**).</br>

##### Inter-run correlations
- [x] Compute inter-run correlation [compute_run_corr.py](analysis_code/preproc/functional/compute_run_corr.py)
- [x] Make inter-run correlations maps with pycortex [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py) or [pycortex_maps_run_corr.sh](analysis_code/preproc/functional/pycortex_maps_run_corr.sh)

##### PRF Gaussian fit
- [x] Create the visual matrix design [vdm_builder.py](analysis_code/postproc/prf/vdm_builder.py)
- [x] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [x] Make pRF maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py) or [pycortex_maps_gridfit.sh](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.sh)

##### PRF ROIs
- [x] Draw ROIs on individual fsnative using Inkscape
- [x] Copy sub-170 containing MMP rois from [RetinoMaps](https://github.com/mszinte/RetinoMaps) project [compute_gauss_gridfit_derivatives.py](https://github.com/mszinte/RetinoMaps/blob/main/analysis_code/atlas/create_170k_mmp_rois_mask.ipynb) and mask areas in the overaly that are not covered by data's field of view.
- [x] Create 170k MMP rois masks [create_mmp_rois_atlas.py](analysis_code/atlas/create_mmp_rois_atlas.py)
- [x] Make ROIS files [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py) or [pycortex_maps_rois.sh](analysis_code/postproc/prf/postfit/pycortex_maps_rois.sh)

##### PRF CSS fit
- [x] CSS fit within the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute CSS statistics [compute_css_stats.py](analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification [css_pcm_sbatch.py](analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make CSS fit derivatives and pcm maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py) or [pycortex_maps_css.sh](analysis_code/postproc/prf/postfit/pycortex_maps_css.sh)
- [x] Make subject WEBGL with pycortex [pycortex_webgl_css.py](analysis_code/postproc/prf/webgl/pycortex_webgl_css.py) or [pycortex_webgl_css.sh](analysis_code/postproc/prf/webgl/pycortex_webgl_css.sh)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](analysis_code/postproc/prf/webgl/publish_webgl.py)
- [x] Make TSV with CSS fit derivatives, pcm and statistics [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make pRF derivatives and pcm main figure TSV [make_rois_fig_tsv.py](analysis_code/postproc/prf/postfit/make_rois_fig_tsv.py)
- [x] Make pRF derivatives and pcm main figures [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py) or [make_rois_fig.sh](analysis_code/postproc/prf/postfit/make_rois_fig.sh)
- [x] Merge all css pycortex and pRF derivatives and pcm main figures [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)

### Group analysis
We run either analysis on the template of the HCP cifti format (**sub-170k**) in which individual results are averaged on a template </br>
or we ran a ROI based group analysis determined individually on subject surfaces fsnative (**group**).</br> 

#### Structural preprocessing
- [x] Compute vertex areas for **sub-170k** [compute_vertex_area.py](analysis_code/preproc/anatomical/compute_vertex_area.py)

#### Functional postprocessing

##### Inter-run correlations
- [x] Compute inter-run correlation for **sub-170k** [compute_run_corr.py](analysis_code/preproc/functional/compute_run_corr.py)
- [x] Make inter-run correlations maps with pycortex for **sub-170k** [pycortex_maps_run_corr.py](analysis_code/preproc/functional/pycortex_maps_run_corr.py)

##### PRF Gaussian fit
- [x] Compute pRF gaussian grid fit derivatives for **sub-170k** [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [x] Make pRF maps with pycortex for **sub-170k**  [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py)

##### PRF ROIs
- [x] Make ROIS files for **sub-170k** [make_rois_img.py](analysis_code/postproc/prf/postfit/make_rois_img.py)
- [x] Create flatmaps of ROIs for **sub-170k** [pycortex_maps_rois.py](analysis_code/postproc/prf/postfit/pycortex_maps_rois.py)

##### PRF CSS fit
- [x] Compute CSS statistics for **sub-170k** [compute_css_stats.py](analysis_code/postproc/prf/postfit/compute_css_stats.py)
- [x] Compute CSS fit derivatives for **sub-170k** [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)
- [x] Compute CSS population cortical magnification for **sub-170k** [css_pcm_sbatch.py](analysis_code/postproc/prf/postfit/css_pcm_sbatch.py)
- [x] Make CSS fit derivatives and pcm maps with pycortex for **sub-170k** [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py)
- [x] Make subject WEBGL with pycortex for **sub-170k** [pycortex_webgl_css.py](analysis_code/postproc/prf/webgl/pycortex_webgl_css.py)
- [x] Edit [index.html](analysis_code/postproc/prf/webgl/index.html) and publish WEBGL on webapp [publish_webgl.py](analysis_code/postproc/prf/webgl/publish_webgl.py)
- [x] Make TSV with CSS fit derivatives, pcm and statistics for **sub-170k** [make_tsv_css.py](analysis_code/postproc/prf/postfit/make_tsv_css.py)
- [x] Make pRF derivatives and pcm main figure TSV for **sub-170k** and **group** [make_roi_fig_tsv.py](analysis_code/postproc/prf/postfit/make_roi_fig_tsv.py)
- [x] Make pRF derivatives and pcm main figure TSV for **sub-170k** and **group** [make_rois_fig.py](analysis_code/postproc/prf/postfit/make_rois_fig.py)
- [x] Merge all css pycortex and pRF derivatives and pcm main figures for **sub-170k** and **group** [merge_fig_css.py](analysis_code/postproc/prf/postfit/merge_fig_css.py)

## On working
### **GLM analysis**
- [x] Run Glm for the differents tasks [glm_sbatch.py](analysis_code/postproc/glm/glm_sbatch.py)
- [x] Compute GLM statistics [compute_glm_stats.py](analysis_code/postproc/glm/postfit/compute_glm_stats.py)

### **Statistical analysis**
- [x] Make final statistique image [stats_final.py](analysis_code/postproc/stats/stats_final.py)
- [x] Average statisticales derivatives from all subjects in 170k template [170k_stats_averaging.py](analysis_code/postproc/stats/170k_stats_averaging.py) 
- [x] Make final satistiques maps maps with pycortex [pycortex_maps_stats_final.py](analysis_code/postproc/stats/pycortex_maps_stats_final.py)