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
</br>

## Data analysis
---

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
- [ ] make cut in the brains for flattening [cortex_cuts.sh](analysis_code/preproc/anatomical/cortex_cuts.sh)
- [ ] flatten the cut brains [flatten_sbatch.py](analysis_code/preproc/anatomical/flatten_sbatch.py)

#### Functional preprocessing
- [ ] fMRIprep [fmriprep_sbatch.py](analysis_code/preproc/functional/fmriprep_sbatch.py)
- [ ] supress bad run [bad_run.py](analysis_code/preproc/functional/bad_run.py)
- [ ] high-pass, z-score, average and leave-one-out average [preproc_end.py](analysis_code/preproc/functional/preproc_end.py)
- [ ] Load freesurfer and execute [pycortex_import.py](analysis_code/preproc/functional/pycortex_import.py): run only [freesurfer_import_pycortex.py](analysis_code/preproc/functional/freesurfer_import_pycortex.py)

### Post-processing

#### PRF analysis
- [ ] create the visual matrix design [vdm_builder.py](analysis_code/postproc/prf/fit/vdm_builder.py)
- [ ] Execute [prf_fit.py](analysis_code/postproc/prf/fit/prf_fit.py) to fit pRF parameters (eccentricity, size, amplitude, baseline, rsquare): run only [submit_fit_jobs.py](analysis_code/postproc/prf/fit/submit_fit_jobs.py)
- [ ] Compute pRF derivatives [compute_derivatives.py](analysis_code/postproc/prf/postfit/compute_derivatives.py)
- [ ] make pycortex maps [pycortex_maps.py](analysis_code/postproc/prf/postfit/pycortex_maps.py)
- [ ] draw ROIs using Inkscape
- [ ] compute population cortical magnification [compute_pcm.py](analysis_code/postproc/prf/postfit/compute_pcm.py)
- [ ] extract ROIs masks [roi_masks.ipynb](analysis_code/postproc/prf/postfit/roi_masks.ipynb) 
- [ ] make pdf files with the maps [pdf_maps.py](analysis_code/postproc/prf/postfit/pdf_maps.py)
- [ ] make webgl with the pycortex dataset [pycortex_maps.py](analysis_code/postproc/prf/webgl/pycortex_webgl.py) 
- [ ] send the files [send_data.sh](analysis_code/postproc/prf/webgl/send_data.sh)

### Main analysis
- [ ] extract all data as pickle files or tsv [make_tsv.ipynb](analysis_code/postproc/prf/postfit/make_tsv.ipynb)
- [ ] think about the individual participants figures
- [ ] Figures and statistics [amblyo_analysis_and_figures.ipynb](analysis_code/postproc/result_analysis/amblyo_analysis_and_figures.ipynb)