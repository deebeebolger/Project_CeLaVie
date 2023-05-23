import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.pyplot as plt

base_dir = '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data'

skool = 'Ecole1'
block = 'RS2'
Sujetnums = ['S18', 'S19', 'S20', 'S21', 'S22']  #

if skool == 'Ecole1':
    data_folder = os.path.join(base_dir, 'E1_PreprocData')
elif skool == 'Ecole2':
    data_folder = os.path.join(base_dir, 'E2_PreprocData')
AllPSD = []

fig, axs = plt.subplots(len(Sujetnums), 1)

for counter, sujs in enumerate(Sujetnums):
    datadir_curr = os.path.join(data_folder, sujs)
    datadir_content = os.listdir(datadir_curr)
    startw = block+'_'+sujs
    currfile_title = [x for x in datadir_content if x.startswith(startw)& x.endswith('.fif')]
    fulldir_curr = os.path.join(datadir_curr, currfile_title[0])
    RawIn = mne.io.read_raw_fif(fulldir_curr, preload=True)
    sfreq = RawIn.info['sfreq']

    RawData = RawIn.get_data()  # Get the rawdata array with structure: electrodes X time
    PSD, freqs = mne.time_frequency.psd_array_welch(RawData[0:127,:], sfreq, fmin=1, fmax=60, n_fft=2048, n_overlap=0, average='mean')
    AllPSD.append(PSD)

    axs[counter].semilogy(freqs, np.transpose(PSD))
    if counter == len(Sujetnums)-1:
        axs[counter].set_xlabel('Frequency (Hz)')
    if counter <len(Sujetnums)-1:
        axs[counter].get_xaxis().set_ticks([])
        axs[counter].spines['bottom'].set_visible(False)
    axs[counter].set_ylabel('log(Power)')
    axs[counter].set_frame_on(1)
    axs[counter].set_title(block + ' ' + sujs)
    axs[counter].spines['top'].set_visible(False)
    axs[counter].spines['right'].set_visible(False)

plt.show()

# Calculate the mean PSD across the participants.
allpsd = np.asarray(AllPSD)
psd_mean = np.mean(allpsd, 0)
plt.semilogy(freqs, np.transpose(psd_mean))
plt.xlabel('Frequency (Hz)')
plt.ylabel('log(power)')
plt.title('Average of '+ skool + ' '+ block + ' ( 5 participants)')
plt.box(1)
plt.show()





