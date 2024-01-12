# RetinoMaps
## About
---
We study cerebral areas implicated in both vision and eye movements..</br>
This repository contain all code allowing us to analyse our dataset [OpenNeuro:DSXXXXX](https://openneuro.org/datasets/dsXXXX).</br>

---
## Authors (alphabetic order): 
---
Uriel LASCOMBES, Guillaume MASSON & Martin SZINTE

### Main dependencies
---
[dcm2niix](https://github.com/rordenlab/dcm2niix); 
[PyDeface](https://github.com/poldracklab/pydeface); 
[fMRIprep](https://fmriprep.org/en/stable/); 
[pRFpy](https://github.com/VU-Cog-Sci/prfpy); 
[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/);
[FFmpeg](https://ffmpeg.org/);
[FSL](https://fsl.fmrib.ox.ac.uk);
[Inkscape](https://inkscape.org/)
[workbench](https://humanconnectome.org/software/connectome-workbench)
</br>


## Data analysis
---

### To Do 
- [ ] Update help section prf_fit
- [ ] In plots.ipynb add boxplot across rois and type of fit
- [ ] Test with prf_fit_test.py grid fit with 100 steps of guassian grid vs iterative no bound vs iterative boudn 5 repeat
- [ ] Code prf_fit_grid: gaussian fit with large grid (100 steps?) and only _avg data (SEE prf_fit_test codes)
- [ ] Plot retinotopy data, and get ROIs with grid fit
- [ ] Code prf_fit_final which takes gaussian grid and run better model (dn or css) only on ROIs vertices

### Pre-processing

#### BIDS
- [x] copy relevant data from PredictEye [copy_data.py](analysis_code/preproc/bids/bids_copy_data.sh) 
- [x] change the 'task' to 'task_condition' coulumn name in event.tsv files to avoid BIDS problems [correct_events_files.ipynb](analysis_code/preproc/bids/correct_events_files.ipynb)
- [x] deface participants t1w image [deface_sbatch.py](analysis_code/preproc/bids/deface_sbatch.py) 
    </br>Note: run script for each subject separately.
- [x] validate bids format [https://bids-standard.github.io/bids-validator/] / alternately, use a docker [https://pypi.org/project/bids-validator/]
    </br>Note: for the webpage, use FireFox and wait for at least 30 min, even if nothing seems to happen.

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
- [x] Load freesurfer and import subject in pycortex db [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)
- [x] High-pass, z-score, run correlations, average and leave-one-out average [preproc_end_sbatch.py](analysis_code/preproc/functional/preproc_end_sbatch.py) 
- [ ] Make timeseries inter-run correlation maps with pycortex [pycortex_corr_maps.py](analysis_code/preproc/functional/pycortex_corr_maps.py)

### Post-processing

#### PRF analysis
- [x] Create the visual matrix design [vdm_builder.ipynb](analysis_code/postproc/prf/fit/vdm_builder.ipynb)
- [ ] Run pRF gaussian grid fit [prf_submit_gridfit_jobs.py](analysis_code/postproc/prf/fit/prf_submit_gridfit_jobs.py)
- [ ] Compute pRF gaussian grid fit derivatives [compute_gauss_gridfit_derivatives.py](analysis_code/postproc/prf/postfit/compute_gauss_gridfit_derivatives.py)
- [ ] Make pRF derivatives maps with pycortex [pycortex_maps_gridfit.py](analysis_code/postproc/prf/postfit/pycortex_maps_gridfit.py)
- [ ] draw ROIs using Inkscape
- [ ] Run pRF CSS fit only on the ROIs [prf_submit_css_jobs.py](analysis_code/postproc/prf/fit/prf_submit_css_jobs.py)
- [ ] Compute pRF CSS fit derivatives [compute_css_derivatives.py](analysis_code/postproc/prf/postfit/compute_css_derivatives.py)


#### GLM analysis
- [ ] Run Glm for the differents tasks [glm_sbatch.py](analysis_code/postproc/glm/glm_sbatch.py)

#### On working
- [ ] compute population cortical magnification [compute_pcm.py](analysis_code/postproc/prf/postfit/compute_pcm.py)
- [ ] extract ROIs masks [roi_masks.ipynb](analysis_code/postproc/prf/postfit/roi_masks.ipynb) 
- [ ] make pdf files with the maps [pdf_maps.py](analysis_code/postproc/prf/postfit/pdf_maps.py)
- [ ] make webgl with the pycortex dataset [pycortex_maps.py](analysis_code/postproc/prf/webgl/pycortex_webgl.py) 
- [ ] send the files [send_data.sh](analysis_code/postproc/prf/webgl/send_data.sh)

### Main analysis
- [ ] extract all data as pickle files or tsv [make_tsv.ipynb](analysis_code/postproc/prf/postfit/make_tsv.ipynb)
- [ ] think about the individual participants figures
- [ ] Figures and statistics [amblyo_analysis_and_figures.ipynb](analysis_code/postproc/result_analysis/amblyo_analysis_and_figures.ipynb)