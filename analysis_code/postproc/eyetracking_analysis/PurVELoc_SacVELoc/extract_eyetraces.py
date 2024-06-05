"""
-----------------------------------------------------------------------------------------
extract_eyetraces.py
-----------------------------------------------------------------------------------------
Goal of the script:
Extract eye traces from edf file and arrange them well for later treatment
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number (sub-01)
sys.argv[2]: task (PurVELoc or SacVELoc)
-----------------------------------------------------------------------------------------
Output(s):
h5 files with loads of data on eye traces across runs
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/RetinoMaps/analysis_code/postproc/eyetracking_analysis
python extract_eyetraces.py sub-02 PurVELoc
-----------------------------------------------------------------------------------------
"""
#%%
# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import platform
import re
import numpy as np
import ipdb
import json
import h5py
import scipy.io
import pandas as pd
import glob 
deb = ipdb.set_trace

from sac_utils import replace_blinks_with_nan
# Get inputs
# ----------
subject = 'sub-25'
task = 'SacVELoc'
#subject = sys.argv[1]
#task = sys.argv[2]

# Define analysis parameters
# --------------------------
with open('behavior_settings.json') as f:
	json_s = f.read()
	analysis_info = json.loads(json_s)

# Get eyelink data
# ----------------

main_dir = analysis_info['main_dir_mac']
edf2asc_dir = analysis_info['edf2asc_dir_mac']
end_file = ''


# Define file list
# ----------------
file_dir = '{exp_dir}/{sub}/ses-02'.format(exp_dir = main_dir, sub = subject)

file_dir_save = '/Users/sinakling/projects/PredictEye/data' #! REPLACE WITH LOCAL DIRECTORY 
add_dir = '/Users/sinakling/projects/PredictEye/locVisEndEMexp/data' #additional files (mat) are stored in extra repo, clone before: https://github.com/mszinte/locVisEndEMexp.git
list_filename = ['{sub}_ses-02_task-{task}_run-01'.format(sub = subject, task = task),
				 '{sub}_ses-02_task-{task}_run-02'.format(sub = subject, task = task)]

data_events = sorted(glob.glob(r'{exp_dir}/{sub}/ses-02/func/*{task}*.tsv'.format(exp_dir=main_dir, sub=subject, task = task)))

assert len(data_events) > 0, "No event files found"
# %%

#%%

# Define experiments details
# --------------------------
if subject[-1]=='t':
    num_run = np.arange(0,analysis_info['num_run_t'],1)
else:
    num_run = np.arange(0,analysis_info['num_run'],1)
num_seq = analysis_info['num_seq']
seq_trs = analysis_info['seq_trs']
eye_mov_seq = analysis_info['eye_mov_seq'] # ?
rads = analysis_info['rads']
pursuits_tr = np.arange(0,seq_trs,2) # ?
saccades_tr = np.arange(1,seq_trs,2) # ?

trials_seq = np.array(analysis_info['trials_seq'])
seq_type   = np.array(analysis_info['seq_type'])


mat_seqDir_filename = '{file_dir}/{sub}/ses-02/add/{sub}_task_dir_sequence.mat'.format(file_dir = add_dir,sub = subject)
matSeqDir = scipy.io.loadmat(mat_seqDir_filename)
matSeqDir = matSeqDir['dir_sequence'][0][0][0].flatten()
matSeqDir = [1 if x==1 else -1 for x in matSeqDir] # x==1 -> ccw, x==2 -> cw


# Extract data
# -------------
dfs = []
record_lines = []

eye_data_runs = []
time_last_run_eye = 0
time_start_eye = np.zeros((1, len(num_run)))
time_end_eye = np.zeros((1, len(num_run)))
time_start_seq = np.zeros((num_seq, len(num_run)))
time_end_seq = np.zeros((num_seq, len(num_run)))
time_start_trial = np.zeros((seq_trs, num_seq, len(num_run)))
time_end_trial = np.zeros((seq_trs, num_seq, len(num_run)))

sampling_rate = 1000

for t_run in np.arange(0, len(num_run), 1):
    print(f"--Extracting data from run: {t_run}--")

    edf_filename_edf = '{file_dir}/func/{filename}_eyeData.edf'.format(file_dir=file_dir, filename=list_filename[t_run])
    edf_filename = '{file_dir}/func/{filename}_eyeData'.format(file_dir=file_dir, filename=list_filename[t_run])

    # Get .msg and .dat file
    if not os.path.exists(edf_filename_edf.replace('.edf', '.msg')):
        os.system('edf2asc{} {} -e -y'.format(end_file, edf_filename_edf))   
        os.rename(edf_filename_edf.replace('.edf', '.asc'), edf_filename_edf.replace('.edf', '.msg'))

    if not os.path.exists(edf_filename_edf.replace('.edf', '.dat')): 
        os.system('edf2asc{end_file} {edf_filename} -s -miss -1.0 -y'.format(end_file=end_file, edf_filename=edf_filename_edf))
        os.rename(edf_filename_edf.replace('.edf', '.asc'), edf_filename_edf.replace('.edf', '.dat'))

    # Initialize lists to store data for each run
    runs = []
    timestamps = []
    sequence_numbers = []
    trial_numbers = []
    onset_offsets = []

    with open('{}.msg'.format(edf_filename)) as msgfid:
        first_last_time = False
        first_time = False
        last_time = False
        seq_num = 0
        trial_num = [0] * num_seq  # List to keep track of trial number for each sequence
        sequences_processed = 0

        while sequences_processed < 9:
            line_read = msgfid.readline()
            if not line_read:
                break  # Exit the loop if the end of file is reached

            la = line_read.split()

            if re.search(r"MSG", line_read):
                if len(la) > 2:
                    if la[2] == 'RECORD_START' and not first_time: 
                        record_lines.append(line_read)
                    if la[2] == 'RECORD_STOP':
                        record_lines.append(line_read)
                if line_read.find('sequence 1 started') != -1 and not first_time:
                    time_start_eye[0, t_run] = float(la[1])
                    first_time = True

                if re.search(r"sequence\s\d+\sstarted", line_read):
                    sequence_match = re.search(r"sequence\s(\d+)\sstarted", line_read)
                    if sequence_match:
                        seq_num = int(sequence_match.group(1)) - 1
                        time_start_seq[seq_num, t_run] = float(la[1])

                if re.search(r"sequence\s\d+\sstopped", line_read):
                    sequence_match_end = re.search(r"sequence\s(\d+)\sstopped", line_read)
                    if sequence_match_end:
                        time_end_seq[seq_num, t_run] = float(la[1])

                if re.search(r"trial\s\d+\sonset", line_read):
                    trial_match = re.search(r"trial\s(\d+)", line_read)
                    trial_num[seq_num] = int(trial_match.group(1)) - 1
                    runs.append(t_run)
                    timestamps.append(float(la[1]))
                    sequence_numbers.append(seq_num)
                    trial_numbers.append(trial_num[seq_num])
                    onset_offsets.append('onset')

                if re.search(r"trial\s\d+\soffset", line_read):
                    trial_match = re.search(r"trial\s(\d+)", line_read)
                    trial_num[seq_num] = int(trial_match.group(1)) - 1
                    runs.append(t_run)
                    timestamps.append(float(la[1]))
                    sequence_numbers.append(seq_num)
                    trial_numbers.append(trial_num[seq_num])
                    onset_offsets.append('offset')
                    trial_num[seq_num] += 1  # Move to the next trial

                    # Update time_start_trial for the next trial if exists
                    next_line = msgfid.readline()
                    if next_line and re.search(r"trial\s\d+\sonset", next_line):
                        trial_match = re.search(r"trial\s(\d+)", next_line)
                        runs.append(t_run)
                        trial_num[seq_num] = int(trial_match.group(1)) - 1
                        timestamps.append(float(next_line.split()[1]))
                        sequence_numbers.append(seq_num)
                        trial_numbers.append(trial_num[seq_num])
                        onset_offsets.append('onset')

            if line_read.find('sequence 9 stopped') != -1 and not last_time:
                time_end_eye[0, t_run] = float(la[1])
                last_time = True
                sequences_processed += 1  # Increment the counter for processed sequences

    # Create DataFrame for the current run
    df = pd.DataFrame({
        "Run": runs,
        'Timestamp': timestamps,
        'Sequence_Number': sequence_numbers,
        'Trial_Number': trial_numbers,
        'Onset_Offset': onset_offsets
    })

    # Append the DataFrame to the list of DataFrames
    dfs.append(df)

    # Load eye coord data
    eye_dat = np.genfromtxt('{}.dat'.format(edf_filename), usecols=(0, 1, 2, 3))  # get pupil size as well 
    eye_data_run = eye_dat[np.logical_and(eye_dat[:, 0] >= time_start_eye[0, t_run], eye_dat[:, 0] <= time_end_eye[0, t_run])]
    
    # Add run number
    eye_data_run = np.concatenate((eye_data_run, np.ones((eye_data_run.shape[0], 1)) * t_run), axis=1)
    
    # col 0 = time
    # col 2 = eye x coord
    # col 3 = eye y coord
    # col 4 = run number

    if t_run == 0:
        eye_data_runs = eye_data_run
        
    else:
        eye_data_runs = np.concatenate((eye_data_runs, eye_data_run), axis=0)

    

# Create final DataFrame outside the loop
final_df = pd.concat(dfs, ignore_index=True)
final_df.to_csv(f'{subject}_{task}_timestamps.csv', index=True)
print("--Done Extracting--")
#%%
final_df = pd.read_csv(f'/Users/sinakling/projects/PredictEye/locVisENDEMexp/stats/behav_analysis/{subject}_{task}_timestamps.csv')
        
#%%
# save information from dataframe in new datastructure
for index, row in final_df.iterrows():
    run = row['Run']
    seq_num = row['Sequence_Number']
    trial_num = row['Trial_Number']
    timestamp = row['Timestamp']
    onset_offset = row['Onset_Offset']
    
    if onset_offset == 'onset':
        time_start_trial[trial_num, seq_num, run] = timestamp
    elif onset_offset == 'offset':
        time_end_trial[trial_num, seq_num, run] = timestamp

print(time_end_seq)

#%%
if task == 'PurVELoc':
    overrides = {
        ('sub-02', (1, 1)): 6234689,
        ('sub-02', (6, 1)): 6369218,
        ('sub-03', (4, 1)): 14429843,
        ('sub-03', (5, 1)): 14468270,
        ('sub-05', (7, 1)): 12244007,
        ('sub-07', (2, 1)): 30505304,
        ('sub-08', (8, 1)): 5077674,
        ('sub-09', (6, 1)): 23904419,
        ('sub-11', (6, 0)): 15150292,
        ('sub-11', (8, 1)): 15810500,
        ('sub-12', (5, 0)): 11959093,
        ('sub-12', (3, 1)): 12491363,
        ('sub-13', (1, 1)): 5670252,
        ('sub-14', (4, 0)): 5313466,
        ('sub-17', (6, 0)): 12099065,
        ('sub-20', (3, 0)): 5159813,
        ('sub-21', (0, 1)): 12303800,
        ('sub-22', (2, 0)): 5086667,
        ('sub-22', (6, 1)): 5805845,
        ('sub-23', (0, 0)): 17512475,
        ('sub-23', (2, 0)): 17570125,
        ('sub-23', (4, 1)): 18223182,
        ('sub-23', (5, 1)): 18261617,
        ('sub-24', (2, 0)): 29411597,
        ('sub-24', (7, 1)): 30175079,
        ('sub-24', (8, 1)): 30194286,
        ('sub-25', (4, 0)): 22254951,
        ('sub-25', (6, 0)): 22312568,
        ('sub-25', (1, 1)): 22779580,
       
       
    }

elif task == "SacVELoc":
     overrides = {
         ('sub-02', (0,0)): 5267590, 
         ('sub-02', (2,0)): 5325230, 
         ('sub-02', (6,1)): 6074810, 
         ('sub-03', (1,1)): 14058272,
         ('sub-05', (6,0)): 11319207,
         ('sub-07', (1,0)): 29577339, 
         ('sub-07', (5,0)): 29692619,
         ('sub-07', (6,1)): 30324008,
         ('sub-08', (8,0)): 4148523, 
         ('sub-08', (3,1)): 4635832,
         ('sub-08', (7,1)): 4751122,
         ('sub-09', (0,0)): 22811811, 
         ('sub-09', (4,0)): 22927201,
         ('sub-09', (7,0)): 23023270,
         ('sub-11', (0,0)): 14607157, 
         ('sub-11', (2,0)): 14664797,
         ('sub-11', (4,1)): 15392820,
         ('sub-14', (6,0)): 5041965,
         ('sub-23', (1,1)): 17848104,
         ('sub-24', (1,0)): 29085278,
         ('sub-25', (7,0)): 22045872,

      
    }


for key, value in overrides.items():
    sub, indices = key
    if sub:
        if sub == subject:
            time_end_seq[indices[0], indices[1]] = value

print(time_end_seq)
#%%
print(time_start_trial)
#%%
for t_run in np.arange(0, len(num_run), 1):
    if task == 'PurVELoc':
            timeStartSeq = time_start_trial[0, seq_type == 1, t_run]

            timeEndSeq = time_end_trial[29, seq_type == 1, t_run]  # last TR with pursuit of the seq


            seq_tmp = []
            for idx in range(len(timeStartSeq)):
                data_idx = (eye_data_run[:, 0] >= timeStartSeq[idx]) & (eye_data_run[:, 0] <= timeEndSeq[idx]) 
                seq = (idx + 1) * np.ones(sum(data_idx)) 
                data = np.c_[eye_data_run[data_idx, :], seq]

                if idx == 0:
                    seq_tmp = data
                else:
                    seq_tmp = np.concatenate((seq_tmp, data), axis=0)

            if t_run == 0:
                eye_data_seqs = seq_tmp 
            else:
                eye_data_seqs = np.concatenate((eye_data_seqs, seq_tmp), axis=0)

            # eye_data_seqs
            # col 0 = time
            # col 1 = eye x coord
            # col 2 = eye y coord
            # col 3 = run number
            # col 4 = seq number

# Put nan for blink time
eye_data_runs_nan_blink = replace_blinks_with_nan(eye_data_runs, sampling_rate)  

if task == 'PurVELoc':
    eye_data_seqs_nan_blink = replace_blinks_with_nan(eye_data_seqs, sampling_rate)
#%%
import matplotlib.pyplot as plt

plt.plot(eye_data_runs_nan_blink[:,1],eye_data_runs_nan_blink[:,2], linestyle = 'dotted')
plt.show()
#%%
# convert to dva
# put eye coordinates in deg from center (flip y axis)
scr_sizeX = 1920
scr_sizeY = 1080
screen_size = np.array([scr_sizeX,scr_sizeY])
ppd = 52.00031911320716



eye_data_runs[:,1] = (eye_data_runs[:,1] - (screen_size[0]/2))/ppd
eye_data_runs[:,2] = -1.0*((eye_data_runs[:,2] - (screen_size[1]/2))/ppd)
eye_data_runs_nan_blink[:,1] = (eye_data_runs_nan_blink[:,1] - (screen_size[0]/2))/ppd
eye_data_runs_nan_blink[:,2] = -1.0*((eye_data_runs_nan_blink[:,2] - (screen_size[1]/2))/ppd)

if task == 'PurVELoc':
    eye_data_seqs[:,1] = (eye_data_seqs[:,1] - (screen_size[0]/2))/ppd
    eye_data_seqs[:,2] = (eye_data_seqs[:,2] - (screen_size[1]/2))/ppd

    eye_data_seqs_nan_blink[:,1] = (eye_data_seqs_nan_blink[:,1] - (screen_size[0]/2))/ppd
    eye_data_seqs_nan_blink[:,2] = (eye_data_seqs_nan_blink[:,2] - (screen_size[1]/2))/ppd

#%%

import matplotlib.pyplot as plt

plt.plot(eye_data_runs_nan_blink[:,1], linestyle = 'dotted')
plt.show()
#%%

# Save all
# --------
session = 'ses-02'
h5_file = "{file_dir}/{sub}_{ses}_task-{task}_eyedata.h5".format(file_dir = file_dir_save, sub = subject, ses = session, task = task)


try: os.system('rm {h5_file}'.format(h5_file = h5_file))
except: pass

h5file = h5py.File(h5_file, "a")

h5file.create_dataset('eye_data_runs',data = eye_data_runs,dtype ='float32')
h5file.create_dataset('eye_data_runs_nan_blink',data = eye_data_runs_nan_blink,dtype ='float32')
h5file.create_dataset('time_start_eye',data = time_start_eye,dtype ='float32')
h5file.create_dataset('time_end_eye',data = time_end_eye,dtype ='float32')
h5file.create_dataset('time_start_seq',data = time_start_seq,dtype ='float32')
h5file.create_dataset('time_end_seq',data = time_end_seq,dtype ='float32')
h5file.create_dataset('time_start_trial',data = time_start_trial,dtype ='float32')
h5file.create_dataset('time_end_trial',data = time_end_trial,dtype ='float32')
h5file.create_dataset('dir_sequence',data = matSeqDir, dtype='int')

if task == 'SacVELoc': 
    h5file.close()


if task == 'PurVELoc':
    h5file.create_dataset('eye_data_seqs',data = eye_data_seqs, dtype='float32')
    h5file.create_dataset('eye_data_seqs_nan_blink',data = eye_data_seqs_nan_blink,dtype ='float32')

    h5file.close()

