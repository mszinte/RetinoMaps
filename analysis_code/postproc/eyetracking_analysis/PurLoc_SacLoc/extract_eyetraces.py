"""
-----------------------------------------------------------------------------------------
extract_eyetraces.py
-----------------------------------------------------------------------------------------
Goal of the script:
Extract eye traces from edf file and arrange them well for later treatment
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number (sub-01)
sys.argv[2]: task (PurLoc or SacLoc)
-----------------------------------------------------------------------------------------
Output(s):
h5 files with loads of data on eye traces across runs
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/RetinoMaps/analysis_code/postproc/eyetracking_analysis
python extract_eyetraces.py sub-02 PurLoc
-----------------------------------------------------------------------------------------
"""
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
import glob 
import math
import pandas as pd
import matplotlib.pyplot as plt

deb = ipdb.set_trace

from sac_utils import replace_blinks_with_nan

# Get inputs
# ----------
subject = sys.argv[1]
task = sys.argv[2]

# Define analysis parameters
# --------------------------
with open('behavior_settings.json') as f:
	json_s = f.read()
	analysis_info = json.loads(json_s)

# Get eyelink data


main_dir = analysis_info['main_dir_mac']
edf2asc_dir = analysis_info['edf2asc_dir_mac']
end_file = ''


# Define file list

file_dir = '{exp_dir}/{sub}/ses-02'.format(exp_dir = main_dir, sub = subject)

file_dir_save = '/Users/sinakling/projects/PredictEye/data' #! REPLACE WITH PERSONAL LOCAL DIRECTORY
list_filename = ['{sub}_ses-02_task-{task}_run-01'.format(sub = subject, task = task),
				 '{sub}_ses-02_task-{task}_run-02'.format(sub = subject, task = task)]

data_events = sorted(glob.glob(r'{exp_dir}/{sub}/ses-02/func/*{task}*.tsv'.format(exp_dir=main_dir, sub=subject, task = task)))

# Define experiments details
# --------------------------
num_run = analysis_info['num_run']
num_seq = analysis_info['num_seq']
seq_trs = analysis_info['seq_trs']
eye_mov_seq = analysis_info['eye_mov_seq']
rads = analysis_info['rads']
pursuits_tr = np.arange(0,seq_trs,2)
saccades_tr = np.arange(1,seq_trs,2)


# Exctract data
# -------------
eye_data_runs = [];
time_last_run_eye   =   0;
time_start_eye = np.zeros((1,num_run))
time_end_eye = np.zeros((1,num_run))
time_start_seq = np.zeros((num_seq,num_run))
time_end_seq = np.zeros((num_seq,num_run))
time_start_trial = np.zeros((seq_trs,num_seq,num_run))
time_end_trial = np.zeros((seq_trs,num_seq,num_run))


# Initialize an empty list to store DataFrames for each run
dfs = []
record_lines = []

for t_run in np.arange(0,num_run,1):
    print(f"--Extracting data from run: {t_run}--")

    edf_filename_edf = '{file_dir}/func/{filename}_eyeData.edf'.format(file_dir=file_dir,filename=list_filename[t_run])
    edf_filename = '{file_dir}/func/{filename}_eyeData'.format(file_dir=file_dir,filename=list_filename[t_run])

    # get .msg and .dat file
    if not os.path.exists(edf_filename_edf.replace('.edf','.msg')):
        os.system('edf2asc{} {} -e -y'.format(end_file, edf_filename_edf))   
        os.rename(edf_filename_edf.replace('.edf','.asc'),edf_filename_edf.replace('.edf','.msg')) 


    if not os.path.exists(edf_filename_edf.replace('.edf','.dat')): 
        os.system('edf2asc{end_file} {edf_filename} -s -miss -1.0 -y'.format(end_file=end_file, edf_filename=edf_filename_edf))
        os.rename(edf_filename_edf.replace('.edf','.asc'),edf_filename_edf.replace('.edf','.dat'))

    # Initialize lists to store data for each run
    runs = []
    timestamps = []
    sequence_numbers = []
    trial_numbers = []
    onset_offsets = []

    msgfid = open('{}.msg'.format(edf_filename))
    first_time = False
    last_time = False
    seq_num = 0
    trial_num = [0] * num_seq  # List to keep track of trial number for each sequence
    sequences_processed = 0

    try:
        while sequences_processed < 9:
            line_read = msgfid.readline()
            if not line_read:
                break  # Exit the loop if the end of file is reached

            la = line_read.split()

            if re.search(r"MSG", line_read):
                if len(la) > 2:
                    if la[2] == 'RECORD_START' and not first_time: 
                        record_lines.append(line_read)
                    if la[2] == 'RECORD_STOP':                        #Eyelink('message', 'RECORD_STOP');
                        record_lines.append(line_read)
                if line_read.find('sequence 1 started') != -1 and not first_time:
                    time_start_eye[0, t_run] = float(la[1])
                    print('first time true')
                    first_time = True

                if re.search(r"sequence\s\d+\sstarted", line_read):
                    sequence_match = re.search(r"sequence\s(\d+)\sstarted", line_read)
                    if sequence_match:
                        seq_num = int(sequence_match.group(1)) - 1
                        time_start_seq[seq_num, t_run] = float(la[1])

                if re.search(r"sequence\s\d+\sstopped", line_read):
                    sequence_match_end = re.search(r"sequence\s(\d+)\sstopped", line_read)
                    if sequence_match_end:
                        #seq_num = int(sequence_match.group(1)) - 1
                        time_end_seq[seq_num, t_run] = float(la[1])

                if re.search(r"trial\s\d+\sonset", line_read):
                    trial_match = re.search(r"trial\s(\d+)", line_read)
                    trial_num[seq_num] = int(trial_match.group(1)) - 1
                    runs.append(t_run)
                    timestamps.append(float(la[1]))
                    sequence_numbers.append(seq_num)
                    trial_numbers.append(trial_num[seq_num])
                    onset_offsets.append('onset')
                    #print('n.e, Trial {} onset at {}'.format(trial_num[seq_num], float(la[1])))

                if re.search(r"trial\s\d+\soffset", line_read):
                    trial_match = re.search(r"trial\s(\d+)", line_read)
                    trial_num[seq_num] = int(trial_match.group(1)) - 1
                    runs.append(t_run)
                    timestamps.append(float(la[1]))
                    sequence_numbers.append(seq_num)
                    trial_numbers.append(trial_num[seq_num])
                    onset_offsets.append('offset')
                    #print('n.e, Trial {} offset at {}'.format(trial_num[seq_num], float(la[1])))
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
                        #print('e., Trial {} onset at {}'.format(trial_num[seq_num], float(next_line.split()[1])))

            if line_read.find('sequence 9 stopped') != -1 and not last_time:
                time_end_eye[0, t_run] = float(la[1])
                print('last time true')
                last_time = True
                sequences_processed += 1  # Increment the counter for processed sequences

    finally:
        msgfid.close()
       

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

	# load eye coord data
    eye_dat = np.genfromtxt('{}.dat'.format(edf_filename),usecols=(0, 1, 2, 3)) # get pupil size as well 
    eye_data_run = eye_dat[np.logical_and(eye_dat[:,0]>=time_start_eye[0,t_run],eye_dat[:,0]<=time_end_eye[0,t_run])]

	# add run number
    eye_data_run = np.concatenate((eye_data_run,np.ones((eye_data_run.shape[0],1))*(t_run)),axis = 1)
	# col 0 = time
	# col 2 = eye x coord
	# col 3 = eye y coord
	# col 4 = run number

    if t_run == 0:
        eye_data_runs = eye_data_run
    else:
        eye_data_runs = np.concatenate((eye_data_runs,eye_data_run), axis=0)

final_df = pd.concat(dfs, ignore_index=True)
#final_df.to_csv(f'{subject}_{task}_timestamps.csv', index = True) 
print("--Done Extracting--")

#final_df = pd.read_csv(f'/Users/sinakling/projects/PredictEye/locEMexp/stats/behav_analysis/{subject}_{task}_timestamps.csv')

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

# override missing sequence times depending on subject 
if task == 'SacLoc':
    overrides = {
        ('sub-02', (0, 1)): 4639938,
        ('sub-03', (2, 0)): 12187356,
        ('sub-04', (0, 0)): 28463815,
        ('sub-04', (8, 0)): 28694385,
        ('sub-04', (2, 1)): 29131423,
        ('sub-04', (4, 1)): 29189069,
        ('sub-05', (1, 0)): 9937947,
        ('sub-05', (7, 1)): 10708999,
        ('sub-06', (7, 0)): 21444493,
        ('sub-07', (2, 1)): 28945732,
        ('sub-08', (0, 1)): 3226434,
        ('sub-09', (7,0)): 21709296,
        ('sub-11', (3,1)): 14063791, 
        ('sub-11', (5,1)): 14121454,
        ('sub-12', (2,0)): 10283144,
        ('sub-13', (0,0)): 3367598, 
        ('sub-13', (6,0)): 3540507,
        ('sub-13', (5,1)): 4169484,
        ('sub-13', (8,1)): 4246314,
        ('sub-14', (2,0)): 3546104,
        ('sub-14', (5,0)): 3642204,
        ('sub-14', (6,0)): 3661433,
        ('sub-14', (1,1)): 4129644,
        ('sub-20', (3,1)): 4217449, 
        ('sub-20', (5,1)): 4275096,
        ('sub-22', (7,1)): 4273942,
        ('sub-23', (6,1)): 16770406,
        ('sub-24', (1,0)): 27812803,
        ('sub-24', (8,0)): 28004889,
        ('sub-24', (0,1)): 28379916,
        ('sub-25', (5,0)): 20755806,



    }

elif task == "PurLoc":
     overrides = {
        ('sub-02', (4, 0)): 4456646,
        ('sub-02', (6, 0)): 4514317,
        ('sub-02', (7, 0)): 4552731,
        ('sub-03', (1, 0)): 12463159,
        ('sub-03', (3, 0)): 12520780,
        ('sub-03', (6, 0)): 12597646,
        ('sub-03', (6, 1)): 13195242,
        ('sub-03', (2, 0)): 12187356,
        ('sub-04', (0, 1)): 28772795,
        ('sub-05', (4, 0)): 10313861,
        ('sub-06', (7, 0)): 21798823,
        ('sub-06', (6, 0)): 21760406,
        ('sub-07', (4, 0)): 28696960,
        ('sub-07', (6, 0)): 28754591,
        ('sub-07', (3, 1)): 29304388,
        ('sub-07', (4, 1)): 29323598,
        ('sub-08', (2, 0)): 2986369,
        ('sub-14', (0, 0)): 3787164,
        ('sub-14', (5, 1)): 4540118,
        ('sub-17', (4, 0)): 10776893,
        ('sub-20', (1, 0)): 3851791,
        ('sub-20', (6, 1)): 4601768,
        ('sub-20', (7, 1)): 4640188,
        ('sub-21', (6, 0)): 10647174,
        ('sub-22', (3, 0)): 3855592,
        ('sub-22', (5, 0)): 3913234,
        ('sub-22', (8, 0)): 3990068,
        ('sub-22', (0, 1)): 4368952,
        ('sub-23', (6, 0)): 16410388,
        ('sub-23', (0, 1)): 16899000,
        ('sub-23', (4, 1)): 17014296,
        ('sub-24', (0, 0)): 28084299,
        ('sub-24', (2, 0)): 28141931,
        ('sub-25', (2, 0)): 20951239,
        ('sub-25', (4, 1)): 21595385,
        ('sub-25', (6, 1)): 21653002,
        ('sub-25', (7, 1)): 21691412,
    }


for key, value in overrides.items():
    sub, indices = key
    if sub:
        if sub == subject:
            time_end_seq[indices[0], indices[1]] = value

print(time_end_seq)

# Extract recording times from msg file 

record_dict = {}
record_start_count = 1
record_stop_count = 1

for line in record_lines:
    parts = line.split()
    timestamp = int(parts[1])
    action = parts[2]

    if action == 'RECORD_START':
        key = f'RECORD_START_{record_start_count}'
        record_dict[key] = timestamp
        record_start_count += 1
    elif action == 'RECORD_STOP':
        key = f'RECORD_STOP_{record_stop_count}'
        record_dict[key] = timestamp
        record_stop_count += 1

print(record_dict)

# Extracting the values from the dictionary and arranging them in the desired structure
record_array = np.array([
    [record_dict['RECORD_START_1'], record_dict['RECORD_STOP_1']],
    [record_dict['RECORD_START_2'], record_dict['RECORD_STOP_2']]
])

# Put nan for blink time
sampling_rate = 1000

eye_data_runs_nan_blink = replace_blinks_with_nan(eye_data_runs,sampling_rate)

#plt.plot(eye_data_runs_nan_blink[:,2], linestyle = 'dotted')
#plt.show()

# Convert to dva and put eye coordinates in deg from center (flip y axis)
# ----------------------------------------------------------------------

scr_sizeX = 1920
scr_sizeY = 1080
screen_size = np.array([scr_sizeX,scr_sizeY])
ppd = 52.00031911320716


eye_data_runs[:,1] = (eye_data_runs[:,1] - (screen_size[0]/2))/ppd
eye_data_runs[:,2] = -1.0*((eye_data_runs[:,2] - (screen_size[1]/2))/ppd)

eye_data_runs_nan_blink[:,1] = (eye_data_runs_nan_blink[:,1] - (screen_size[0]/2))/ppd
eye_data_runs_nan_blink[:,2] = -1.0*((eye_data_runs_nan_blink[:,2] - (screen_size[1]/2))/ppd)

# get amplitude sequence from event files
dfs = []
legend_amp = {1: 4, 2: 6, 3: 8, 4: 10, 5: "none"}
    
for file_path in data_events:
	df = pd.read_csv(file_path, sep='\t')
	dfs.append(df)

appended_df = pd.concat(dfs, ignore_index=True)
amp_sequence_ev = list(appended_df['eyemov_amplitude'])

amp_sequence = [legend_amp[val] if not math.isnan(val) else float('nan') for val in amp_sequence_ev]

# Save all
# --------
#%%
h5_file = '{file_dir}/{sub}_task-{task}_eyedata.h5'.format(file_dir = file_dir_save, sub = subject, task = task)


try:
    os.system(f'rm "{h5_file}"')  
except:
    pass

h5file = h5py.File(h5_file, "a")


# Creating datasets
h5file.create_dataset(f'eye_data_runs', data=eye_data_runs, dtype='float32')
h5file.create_dataset(f'eye_data_runs_nan_blink', data=eye_data_runs_nan_blink, dtype='float32')
h5file.create_dataset(f'time_start_eye', data=time_start_eye, dtype='float32')
h5file.create_dataset(f'time_end_eye', data=time_end_eye, dtype='float32')
h5file.create_dataset(f'time_start_seq', data=time_start_seq, dtype='float32')
h5file.create_dataset(f'time_end_seq', data=time_end_seq, dtype='float64')
h5file.create_dataset(f'time_start_trial', data=time_start_trial, dtype='float32')
h5file.create_dataset(f'time_end_trial', data=time_end_trial, dtype='float32')
h5file.create_dataset(f'amp_sequence', data=amp_sequence_ev, dtype='float32')
h5file.create_dataset(f'recording_times', data=record_array, dtype='float32')


h5file.close()


