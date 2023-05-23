import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.pyplot as plt

data_folder_E1 = '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data/E1_PreprocData'
data_folder_E2 = '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data/E2_PreprocData'

skool = 'Ecole1'
blocknom = 'RS1'
Sujetnum = 'S18'

if skool == 'Ecole1':
    fulldir_curr = os.path.join(data_folder_E1, Sujetnum)
    fulldir_content = os.listdir(fulldir_curr)
elif skool == 'Ecole2':
    fulldir_curr = os.path.join(data_folder_E2, Sujetnum)
    fulldir_content = os.listdir(fulldir_curr)


filenom = [x for x in fulldir_content if x.startswith(blocknom) & x.endswith('.fif')]
fullpath_data = os.path.join(fulldir_curr, filenom[0])
RawIn = mne.io.read_raw_fif(fullpath_data, preload=False)

events_curr = mne.find_events(RawIn, initial_event=True, stim_channel=None)

# Visualize the continuous data
RawIn.plot_sensors()
RawIn.plot(block=True, title= 'Continuous Data : ' + filenom[0])

# Resample the data to 128Hz (for microstate analysis).
srate_div = 4
srate = RawIn.info['sfreq']
srate_new = srate/srate_div
RawIn_rs = RawIn.copy().resample(sfreq=srate_new)

RawInLP_pick = RawIn_lpass.copy().pick_types(meg=False, eeg=True, stim=False, eog=False, misc=False, exclude=[])
RawIn_pick = RawIn.copy().pick_types(meg=False, eeg=True, stim=False, eog=False, misc=False, exclude=[])
RawInLP_pick.plot(block=True, title= 'Continuous Data low-pass filtered (40Hz) : ' + filenom[0], n_channels=40, bad_color='r', butterfly=False)

fig, axs = plt.subplots(2, 1)
for scnt, axs_curr in enumerate(axs):
 mne.viz.plot_raw_psd(RawInLP_pick, fmin=0.1, fmax=60, ax=axs_curr)




