import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.pyplot as plt

def create_dir(save_dir, sujcurr):

    # If the folder for the current subjects of interest do not exist in the save directory, create a folder for current subject.
    curr_savepath = os.path.join(save_dir, sujcurr)

    if os.path.exists(curr_savepath):
        print(f'The save directory {curr_savepath} already exists.\n')
    else:
        print(f'The save directory {curr_savepath} does not exist. Creating...\n')
        os.makedirs(curr_savepath)

    return curr_savepath

def create_filter(f_p, sfreq, dataIn):
    """
    Function to apply a long-duration, FIR low-pass filter with a steep cutoff.
    :param f_p: low-pass cutoff.
    :param sfreq: sampling frequency.
    :param dataIn: input data.
    :return: datafilt - filtered data.
    """
    transBW = 0.5
    f_s = f_p + transBW
    filtdur = 10
    filtlen = int(sfreq*filtdur)

    datafilt = dataIn.copy().filter(None, f_p, h_trans_bandwidth=transBW, filter_length= '%ss' % filtdur)
    # datafilt = mne.filter.filter_data(dataIn, sfreq, l_freq=None, h_freq=f_p, picks=None, filter_length= '%ss' % filtdur,
    #                        l_trans_bandwidth=transBW)

    return datafilt


base_dir = '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data'
Data_folder_E1= '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data/Ecole1'
Data_folder_E2= '/Users/bolger/Documents/work/Projects/Project_CeLaVie/Data/Ecole2'

skool = 'Ecole2'
blocknom = 'RS1'
Sujetnums = ['S42']  #,'S19','S20','S21', 'S22'
Data_folder_curr = os.path.join(base_dir, skool)

if skool == 'Ecole1':
    Save_folder = os.path.join(base_dir, 'E1_PreprocData')
elif skool == 'Ecole2':
    Save_folder = os.path.join(base_dir, 'E2_PreprocData')

eog_channels = ['EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8']
misc_channels = ['GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp']
srate_div = 32  # To have sampling frequency of 128Hz (required for micro-state analysis)

for cntr, suj in enumerate(Sujetnums):

        Savefolder_curr = create_dir(Save_folder, suj)  # Call of function to
        dir_suj = os.path.join(Data_folder_curr, suj)
        dirsuj_content = os.listdir(dir_suj)
        fnom_curr = blocknom+'_'+suj[1:]+'.bdf'
        fulldir   = os.path.join(dir_suj, fnom_curr)
        RawIn     = mne.io.read_raw_bdf(fulldir, eog=eog_channels , misc=misc_channels, stim_channel='auto', preload=True)

        # Print current data information.
        channoms = RawIn.info['ch_names']
        print(channoms)
        srate = RawIn.info['sfreq']
        print(f'The sampling rate of the current dataset is {srate}Hz')

        # Apply Montage and plot sensor layout.
        montage = mne.channels.make_standard_montage('biosemi128')
        RawIn.set_montage(montage)
        RawIn.plot_sensors(show_names=True)

        # Carry out basic pre-processing steps.
        RawIn_ref = mne.set_eeg_reference(RawIn, ref_channels='average')  # Apply the average reference
        new_srate = srate / srate_div
        print(f'New sampling rate is {new_srate}Hz')
        RawIn_rs = RawIn_ref[0].copy().resample(sfreq=new_srate)
        RawIn_hpass = RawIn_rs.copy().filter(l_freq=0.1, h_freq=None)   # Highpass the data.
        #RawIn_lpass = RawIn_hpass.copy().filter(l_freq=None, h_freq=40, ) # Lowpass the data.
        RawIn_lpass = create_filter(40, new_srate, RawIn_hpass)  # Call of filter to apply a FIR filter with steep transition BW.

        # Visualize the continuous data to mark bad electrodes.
        RawIn_lpass.plot(block=True, n_channels=40, bad_color='r', butterfly=False)  # Plot the continuous data.

        # Save the preprocessed object as .fif.
        save_nom = blocknom+'_'+suj+'-ref-rs-hp-lp_v2.fif'
        save_fullfile = os.path.join(Savefolder_curr, save_nom)
        RawIn_lpass.save(save_fullfile, overwrite=True)

        # Extract the event information from the stim channel.
        events_curr = mne.find_events(RawIn_lpass, initial_event=True, stim_channel=None)
        fsplit = save_nom.split('-')
        if len(events_curr)>1:

            event_fnom = fsplit[0] + '-events.txt'
            save_events = os.path.join(save_fullfile, event_fnom)
            mne.write_events(save_events, events_curr, overwrite=True)

        AllPicks = []
        RawInLP_pick = RawIn_lpass.copy().pick_types(meg=False, eeg=True, stim=False, eog=False, misc=False, exclude=[], selection=channoms[0:128])
        AllPicks.append(RawInLP_pick)
        RawInLP_pick2 = RawIn_lpass.copy().pick_types(meg=False, eeg=True, stim=False, eog=False, misc=False, exclude='bads')
        AllPicks.append(RawInLP_pick2)

        fig, axs = plt.subplots(2, 1)
        for scnt, axs_curr in enumerate(axs):
            print(scnt)
            L = len(AllPicks[scnt].info['ch_names'])
            if L == 128:
                picks = channoms[0:128]
                currtitle = fsplit[0]+' : '+str(L)+' channels'
            else:
                picks = 'all'
                currtitle = fsplit[0] + ' : ' + str(L) + ' channels'
            mne.viz.plot_raw_psd(AllPicks[scnt], fmin=0.1, fmax=60, picks=picks, exclude=[],ax=axs_curr)
            axs_curr.set_title(currtitle)



