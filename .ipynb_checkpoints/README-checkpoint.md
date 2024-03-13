# RetinoMaps
## About
---
We study cerebral areas implicated in both vision and eye movements..</br>
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
- [x] Make timeseries inter-run correlation maps with pycortex [pycortex_corr_maps.py](analysis_code/preproc/functional/pycortex_corr_maps.py)

## Post-processing
---
### **pRF analysis**
- [x] Create the visual matrix design [vdm_builder.ipynb](analysis_code/postproc/prf/fit/vdm_builder.ipynb)

#### Gaussian fit
- [x] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [x] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py) 
- [x] Make pRF derivatives maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py)

#### Rois
- [x] Draw ROIs on fsnative using Inkscape
- [x] Import ROIs from the MMP atlas for the hcp 170k template [creat_hcp_rois.ipynb](analysis_code/atlas/creat_hcp_rois.ipynb)
- [x] Creat 170k MMP rois mask [creat_170k_mmp_rois_mask.ipynb](analysis_code/atlas/creat_170k_mmp_rois_mask.ipynb)

#### Css fit
- [x] Run pRF CSS fit only on the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [x] Compute pRF CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)

### pCM analysis
- [x] Compute population cortical magnification [pcm_sbatch.py](analysis_code/postproc/pcm/pcm_sbatch.py)

### Final pycortex maps 
- [x] Make pRF derivatives and pcm maps with pycortex [pycortex_maps_css.py](analysis_code/postproc/prf/postfit/pycortex_maps_css.py)

#### Final ploting
- [x] Make sub-all tsv [make_sub_all_tsv.py](analysis_code/postproc/prf/postfit/make_sub_all_tsv.py)
- [x] Make figures for all subjects [finals_figures.py](analysis_code/postproc/prf/postfit/finals_figures.py)

### **GLM analysis**
- [x] Run Glm for the differents tasks [glm_fit.py](analysis_code/postproc/glm/glm_fit.py)
- [x] Compute Glm derivatives [compute_glm_derivatives.py](analysis_code/postproc/glm/compute_glm_derivatives.py)
- [x] Make GLM derivatives maps with pycortex [pycortex_maps_glm_final.py](analysis_code/postproc/glm/pycortex_maps_glm_final.py)

### **Sub all analysis**
- [x] Merge all subjects pycortex maps [pdf_maps.py](analysis_code/postproc/pdf_maps.py)

## On working
---
### pCM analysis
- [ ] Compute population cortical magnification [compute_pcm.py](analysis_code/postproc/pcm/compute_pcm.py)
 
### webgl
- [ ] Make webgl with the pycortex dataset [pycortex_maps.py](analysis_code/postproc/prf/webgl/pycortex_webgl.py)
- [ ] Send the files [send_data.sh](analysis_code/postproc/prf/webgl/send_data.sh)

