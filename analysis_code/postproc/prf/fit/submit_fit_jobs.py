"""
-----------------------------------------------------------------------------------------
submit_fit_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create jobscript to fit pRF
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
>> cd to function
>> python fit/submit_fit_jobs.py [pp directory] [subject]
-----------------------------------------------------------------------------------------
Exemple:
1. cd to function
>> cd ~/projects/amblyo_prf/analysis_code/postproc/prf/fit
2. run python command
python submit_fit_jobs.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
python submit_fit_jobs.py /scratch/mszinte/data amblyo_prf sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import numpy as np
import os
import json
import sys
import glob
import nibabel as nb
import datetime
import ipdb
deb = ipdb.set_trace

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Cluster settings
fit_per_hour = 15000.0
nb_procs = 8

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
prf_dir = "{}/{}/prf".format(pp_dir, subject)
os.makedirs(prf_dir, exist_ok=True)
prf_fit_dir = "{}/{}/prf/fit".format(pp_dir, subject)
os.makedirs(prf_fit_dir, exist_ok=True)
prf_jobs_dir = "{}/{}/prf/jobs".format(pp_dir, subject)
os.makedirs(prf_jobs_dir, exist_ok=True)
prf_logs_dir = "{}/{}/prf/log_outputs".format(pp_dir, subject)
os.makedirs(prf_logs_dir, exist_ok=True)

# define permission cmd
chmod_cmd = "\nchmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "\nchgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# Define fns (filenames)
vdm_fn = "{}/{}/derivatives/vdm/vdm.npy".format(main_dir, project_dir)
pp_avg_fns = glob.glob("{}/{}/func/fmriprep_dct_avg/*avg*.nii.gz".format(pp_dir,subject))
for fit_num, pp_avg_fn in enumerate(pp_avg_fns):
    
    input_fn = pp_avg_fn
    fit_fn = "{}/{}_prf-fit.nii.gz".format(prf_fit_dir, os.path.basename(pp_avg_fn)[:-7])
    pred_fn = "{}/{}_prf-pred.nii.gz".format(prf_fit_dir, os.path.basename(pp_avg_fn)[:-7])

    if os.path.isfile(fit_fn):
        if os.path.getsize(fit_fn) != 0:
            print("output file {} already exists: aborting analysis".format(fit_fn))
            exit()

    data = nb.load(input_fn).get_fdata()
    data_var = np.var(data,axis=-1)
    mask = data_var!=0.0    
    num_vox = mask[...].sum()
    job_dur_obj = datetime.timedelta(hours=np.ceil(num_vox/fit_per_hour))
    job_dur = "{:1d}-{:02d}:00:00".format(job_dur_obj.days,divmod(job_dur_obj.seconds,3600)[0])

    # create job shell
    slurm_cmd = """\
#!/bin/bash
#SBATCH -p skylake
#SBATCH -A b327
#SBATCH --nodes=1
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={job_dur}
#SBATCH -e {log_dir}/{sub}_fit_{fit_num}_%N_%j_%a.err
#SBATCH -o {log_dir}/{sub}_fit_{fit_num}_%N_%j_%a.out
#SBATCH -J {sub}_fit_{fit_num}\n\n""".format(
    nb_procs=nb_procs, log_dir=prf_logs_dir, job_dur=job_dur, sub=subject, fit_num=fit_num)

    # define fit cmd
    fit_cmd = "python prf_fit.py {} {} {} {} {} {}".format(
        subject, input_fn, vdm_fn, fit_fn, pred_fn, nb_procs)

    # create sh
    sh_fn = "{}/jobs/{}_prf_fit-{}.sh".format(prf_dir,subject,fit_num)

    of = open(sh_fn, 'w')
    of.write("{}{}{}{}".format(slurm_cmd, fit_cmd, chmod_cmd, chgrp_cmd))
    of.close()

    # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))