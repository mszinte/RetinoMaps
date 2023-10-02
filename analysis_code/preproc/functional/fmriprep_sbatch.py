"""
-----------------------------------------------------------------------------------------
fmriprep_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run fMRIprep on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: server nb of hour to request (e.g 10)
sys.argv[5]: anat only (1) or not (0)
sys.argv[6]: use of aroma (1) or not (0)
sys.argv[7]: use Use fieldmap-free distortion correction
sys.argv[8]: skip BIDS validation (1) or not (0)
sys.argv[9]: save cifti hcp format data with 170k vertices
sys.argv[10]: dof number (e.g. 12)
sys.argv[11]: email account
sys.argv[12]: data group (e.g. 327)
sys.argv[13]: project name (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
preprocessed files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/stereo_prf/analysis_code/preproc/functional
2. run python command
python fmriprep_sbatch.py [main directory] [project name] [subject num]
                          [hour proc.] [anat only] [aroma] [fmapfree] 
                          [skip bids validation] [cifti] [dof] [email account] [group]
-----------------------------------------------------------------------------------------
Exemple:
python fmriprep_sbatch.py /scratch/mszinte/data amblyo_prf sub-01 15 1 0 1 0 1 12 
                            martin.szinte@univ-amu.fr 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# imports modules
import sys
import os
import json
import ipdb
opj = os.path.join
deb = ipdb.set_trace

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[-2:]
hour_proc = int(sys.argv[4])
anat = int(sys.argv[5])
aroma = int(sys.argv[6])
fmapfree = int(sys.argv[7])
skip_bids_val = int(sys.argv[8])
hcp_cifti_val = int(sys.argv[9])
dof = int(sys.argv[10])
email = sys.argv[11]
group = sys.argv[12]
project_name = sys.argv[13]

# Define cluster/server specific parameters
cluster_name  = 'skylake'
singularity_dir = "{main_dir}/{project_dir}/code/singularity/fmriprep-20.2.3.simg".format(
    main_dir=main_dir, project_dir=project_dir)
nb_procs = 32
memory_val = 100
log_dir = "{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs".format(
    main_dir=main_dir, project_dir=project_dir)

# special input
anat_only, use_aroma, use_fmapfree, anat_only_end, use_skip_bids_val, \
    hcp_cifti, tf_export, tf_bind = '','','','','', '', '', ''
if anat == 1:
    anat_only = ' --anat-only'
    anat_only_end = '_anat'
    nb_procs = 8

if aroma == 1:
    use_aroma = ' --use-aroma'

if fmapfree == 1:
    use_fmapfree= ' --use-syn-sdc'

if skip_bids_val == 1:
    use_skip_bids_val = ' --skip_bids_validation'

if hcp_cifti_val == 1:
    tf_export = 'export SINGULARITYENV_TEMPLATEFLOW_HOME=/opt/templateflow'
    tf_bind = "{main_dir}/{project_dir}/code/singularity/fmriprep_tf/:/opt/templateflow".format(
        main_dir=main_dir, project_dir=project_dir)
    hcp_cifti = ' --cifti-output 170k'
    
# define SLURM cmd
slurm_cmd = """\
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -p {cluster_name}
#SBATCH --mail-user={email}
#SBATCH -A {proj_name}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.out
#SBATCH -J {subject}_fmriprep{anat_only_end}
#SBATCH --mail-type=BEGIN,END\n\n{tf_export}
""".format(proj_name=proj_name, nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, anat_only_end=anat_only_end, memory_val=memory_val,
           log_dir=log_dir, email=email, tf_export=tf_export,
           cluster_name=cluster_name)

# define singularity cmd
singularity_cmd = "singularity run --cleanenv -B {tf_bind} -B {main_dir}:/work_dir {simg} --fs-license-file /work_dir/{project_dir}/code/freesurfer/license.txt /work_dir/{project_dir}/ /work_dir/{project_dir}/derivatives/fmriprep/ participant --participant-label {sub_num} --bold2t1w-dof {dof} --bold2t1w-init header --output-spaces T1w fsaverage{hcp_cifti} --low-mem --mem-mb {memory_val}000 --nthreads {nb_procs:.0f} {anat_only}{use_aroma}{use_fmapfree}{use_skip_bids_val}".format(
        tf_bind=tf_bind, main_dir=main_dir, project_dir=project_dir,
        simg=singularity_dir, sub_num=sub_num, nb_procs=nb_procs,
        anat_only=anat_only, use_aroma=use_aroma, use_fmapfree=use_fmapfree,
        use_skip_bids_val=use_skip_bids_val, hcp_cifti=hcp_cifti, memory_val=memory_val,
        dof=dof)

# define permission cmd
chmod_cmd = "\nchmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "\nchgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# create sh folder and file
sh_fn = "{main_dir}/{project_dir}/derivatives/fmriprep/jobs/sub-{sub_num}_fmriprep{anat_only_end}.sh".format(
        main_dir=main_dir, sub_num=sub_num,
        project_dir=project_dir, anat_only_end=anat_only_end)


os.makedirs("{main_dir}/{project_dir}/derivatives/fmriprep/jobs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)
os.makedirs("{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs".format(
                main_dir=main_dir,project_dir=project_dir), exist_ok=True)

of = open(sh_fn, 'w')
of.write("{slurm_cmd}{singularity_cmd}{chmod_cmd}{chgrp_cmd}".format(
    slurm_cmd=slurm_cmd, singularity_cmd=singularity_cmd, 
    chmod_cmd=chmod_cmd, chgrp_cmd=chgrp_cmd))
of.close()

# Submit jobs
print("Submitting {sh_fn} to queue".format(sh_fn=sh_fn))
os.chdir(log_dir)
os.system("sbatch {sh_fn}".format(sh_fn=sh_fn))