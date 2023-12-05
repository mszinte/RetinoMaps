"""
-----------------------------------------------------------------------------------------
glm_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run computation of population cortical magnification
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/glm
2. run python command
>> python glm_sbatch.py [main directory] [project name] [subject num] [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
python glm_sbatch.py /scratch/mszinte/data RetinoMaps sub-02 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import json
import os
import sys
import ipdb
deb = ipdb.set_trace

# define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]

# Define cluster/server specific parameters
cluster_name  = analysis_info['cluster_name']
proj_name = analysis_info['project_name']
nb_procs = 8
memory_val = 48
hour_proc = 5

# set folders
log_dir = "{}/{}/derivatives/pp_data/{}/log_outputs".format(main_dir, project_dir, subject)
job_dir = "{}/{}/derivatives/pp_data/{}/jobs".format(main_dir, project_dir, subject)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_glm_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_glm_%N_%j_%a.out
#SBATCH -J {subject}_glm
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=log_dir)
    
glm_fit_cmd = "python glm_fit.py {} {} {} {}".format(main_dir, project_dir, subject, group)


# Define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# create sh fn
sh_fn = "{}/{}_glm.sh".format(job_dir, subject)



of = open(sh_fn, 'w')
of.write("{} \n{} \n{} \n{}".format(slurm_cmd,glm_fit_cmd,chmod_cmd,chgrp_cmd))
of.close()


# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))