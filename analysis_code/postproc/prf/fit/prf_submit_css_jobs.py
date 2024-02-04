"""
-----------------------------------------------------------------------------------------
prf_submit_css_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create jobscript to make a css fit for pRF
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit
2. run python command
python prf_submit_css_jobs.py [main directory] [project name] [subject num] [group] [server project] [memory]
-----------------------------------------------------------------------------------------
Exemple:
python prf_submit_css_jobs.py /scratch/mszinte/data RetinoMaps sub-09 327 b327 100
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports

import os
import json
import sys
import glob
import ipdb

deb = ipdb.set_trace

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = sys.argv[6]
hour_proc = 6

# Cluster settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
cluster_name  = analysis_info['cluster_name']
nb_procs = analysis_info['nb_procs_fit_prf_css']
fit_per_hour = analysis_info['fit_per_hour_prf'] 

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)


# define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)
wb_command_cmd = 'export PATH=$PATH:/scratch/mszinte/data/RetinoMaps/code/workbench/bin_rh_linux64'

# Define fns (filenames)
dct_avg_gii_fns = "{}/{}/fsnative/func/fmriprep_dct_loo_avg/*_task-pRF_*avg*.func.gii".format(pp_dir,subject)


pp_fns=  glob.glob(dct_avg_gii_fns) 



for fit_num, pp_fn in enumerate(pp_fns):
    prf_dir = "{}/{}/fsnative/prf".format(pp_dir, subject)
    os.makedirs(prf_dir, exist_ok=True)

    prf_jobs_dir = "{}/{}/fsnative/prf/jobs".format(pp_dir, subject)
    os.makedirs(prf_jobs_dir, exist_ok=True)
    prf_logs_dir = "{}/{}/fsnative/prf/log_outputs".format(pp_dir, subject)
    os.makedirs(prf_logs_dir, exist_ok=True)
    
    
    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_css_fit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_css_fit_%N_%j_%a.out
#SBATCH -J {subject}_css_fit
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=prf_logs_dir)

    # define fit cmd
    fit_cmd = "python prf_cssfit.py {} {} {} {} {} ".format(
        main_dir, project_dir, subject, pp_fn, nb_procs )
    
    # create sh
    sh_fn = "{}/jobs/{}_prf_css_fit-{}.sh".format(prf_dir,subject,fit_num)

    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{} \n{}".format(slurm_cmd, wb_command_cmd, fit_cmd, chmod_cmd, chgrp_cmd))
    of.close()

    #Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))
    

    


    






