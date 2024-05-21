%% Test script to implement the RELAX pipeline.
%%***************************************************
% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):
% use fileseparators 'filesep' for increased compatability if necessary (replace the \ with a filesep [outside of quotes])

% MATLAB signal processing toolbox (from MATLAB website)
% MATLAB statistics and machine learning toolbox (from MATLAB website)

% PREP pipeline to reject bad electrodes (install plugin to EEGLAB, or via the github into the EEGLAB plugin folder):
% Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for large-scale EEG analysis. Frontiers in neuroinformatics, 9, 16.
% http://vislab.github.io/EEG-Clean-Tools/ or install via EEGLAB extensions

% Specify the MWF path:
% https://github.com/exporl/mwf-artifact-removal
% Somers, B., Francart, T., & Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of neural engineering, 15(3), 036007.

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.

% fastica:
% http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml

% ICLabel in your eeglab folder as a plugin or via the github:
% https://github.com/sccn/ICLabel
function [RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants,...
    RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
    RELAXProcessing_ICA_AllParticipants, RELAXProcessingRoundThreeAllParticipants,...
    RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = CLV_RunRELAX()

%% Check dependencies are installed.
PrepFolderLocation = '/Users/bolger/Documents/MATLAB/PrepPipeline/';
addpath(genpath(PrepFolderLocation), '-end'); % Use genpath to add with subfolders.




%% Call of function to create the RELAX configuration structure.

RELAX_cfg = CLV_CreateRELAX();

%% Check that the folder to accept the processed data has been exists, if not, create it.
%  The mkdir function appears not to function on m3 mac.

if ~exist(RELAX_cfg.OutputPath, 'dir')
    fprintf("Creating folder entitled: RELAXProcessed...\n")
    [status, msg] = mkdir(RELAX_cfg.OutputPath);
else
    fprintf("The folder RELAXProcessed already exists. No need to create a directory. \n")
end

%% Call of function to carry out different processing steps, based on RELAX_Wrapper().
%  Find the number and type of files in the rawdata folder
%  (RELAX_cfg.myPathRaw)

[filename, filepath] = uigetfile({'*.bdf; *.set' 'BDF or eeglab set file';'*.bdf' 'BDF'; '*.set' 'eeglab set file'},...
    'Choose a BDF or .set data file ', 'multiselect', 'on');
RELAX_cfg.filename = filename;
RELAX_cfg.filepath = filepath;

RELAX_cfg.FileNumber = numel(string(filename));
if numel(string(filename))==1
    RELAX_cfg.SingleFile = 1;
    RELAX_cfg.FilesToProcess = 1;
    filename = {filename};
else
    RELAX_cfg.SingleFile = 0;                           % In the case of multiple files.
    RELAX_cfg.FilesToProcess = numel(string(filename)); % Number of datasets to process.
end

%% Loop through the selected datasets to carry out preprocessing.

for fcounter = 1:RELAX_cfg.FilesToProcess

    fprintf("Preprocessing dataset number %d and title %s.\n ", fcounter, filename{1,fcounter});
    [~,FileName,ext] = fileparts(filename{1, fcounter});

    %% Load in the dataset corresponding to the current filename. It is
    % possible to load in either *.bdf or .set format.
    
    if strcmp(ext, '.set')
        EEG = pop_load(filename{1, fcounter}, filepath);
    elseif strcmp(ext, '.bdf')
        [EEG, ~,~] = pop_biosig(fullfile(filepath,filename{1,fcounter}));
    end

    %% Add RELAX processing information to current EEG structure.

    EEG.RELAXProcessing.aFileName=cellstr(FileName);
    EEG.RELAXProcessingExtremeRejections.aFileName=cellstr(FileName);
    EEG.setname = FileName;
    EEG.filepath = filepath;
    EEG.RELAX.Data_has_been_averagerereferenced=0;
    EEG.RELAX.Data_has_been_cleaned=0;
    RELAX_cfg.ms_per_sample=(1000/EEG.srate);

    savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_cfg'];
    save(savefileone,'RELAX_cfg');

     if RELAX_cfg.ms_per_sample<0.7
        warning('The sampling rate for this file is quite high. Depending on your processing power, RELAX may run slowly or even stall. RELAX was validated using 1000Hz sampling rates.');
        warning('To address this, you could downsample your data with: EEG = pop_resample( EEG, 1000), then save the downsampled data prior to running RELAX');
     end

     %% Add the channel locations to the EEG.chanlocs field of the current dataset.
    
     load(RELAX_cfg.caploc,"chanLocations128")
     for locnt = 1:length(chanLocations128)
         EEG.chanlocs(locnt).labels = chanLocations128(locnt).labels;
         EEG.chanlocs(locnt).theta = chanLocations128(locnt).theta;
         EEG.chanlocs(locnt).radius = chanLocations128.radius;
         [EEG.chanlocs(locnt).X, EEG.chanlocs(locnt).Y, EEG.chanlocs(locnt).Z] = deal(chanLocations128(locnt).X, chanLocations128(locnt).Y,...
             chanLocations128(locnt).Z);
         [EEG.chanlocs(locnt).sph_theta, EEG.chanlocs(locnt).sph_phi, EEG.chanlocs(locnt).sph_radius] = deal(chanLocations128(locnt).sph_theta, ...
             chanLocations128(locnt).sph_phi, chanLocations128(locnt).sph_radius);
         EEG.chanlocs(locnt).urchan = chanLocations128(locnt).urchan;   
     end
    % Plot topography showing les channel layout.
    fprintf('Plotting topography of the currently applied channel layout.\n')
    figure; topoplot([], chanLocations128, 'electrodes', 'labels')
    % Save a *.set file with the channel locations added.
    fprintf('Saving current dataset, %s, with channel locations added as *.set file in %s...\n', FileName, RELAX_cfg.OutputPath);
    EEG = pop_saveset( EEG, 'filename',FileName, 'filepath',RELAX_cfg.OutputPath);

    %% Delete those channels marked for deletion or considered irrelevant for the current study.
    %  Should we delete the auxiliary channels (channels 137 tp 143). 
    
    if ~isempty(RELAX_cfg.ElectrodesToDelete)
        EEG=pop_select(EEG,'nochannel',RELAX_cfg.ElectrodesToDelete);
        EEG = eeg_checkset( EEG );
    else
        fprintf('No electrodes marked for rejection at this point.\n')
    end
    EEG.allchan=EEG.chanlocs; % take list of all included channels before any rejections
    
    % Save *.set file with EEG channels removed.
    FileName_nochan = [FileName, '-nochan'];
    EEG = pop_saveset(EEG, 'filename', FileName_nochan, 'filepath', RELAX_cfg.OutputPath);

    %% Band Pass filter the continuous data 
    %  Need to verify the high-pass and low-pass filter limits applied when
    %  analysing micro-states.
    %  The RELAX toolbox applies a 4th order, zero-phase Butterworth filter. Note that
    %  this type of filter should, ideally, not be used prior to ERP
    %  analysis especially with a highpass cutoff > 0.2Hz, this can lead to
    %  distortions. 
    %  1. A notch filter is applied to reduce line noise (47-53Hz).
    %  2  A bandpass filter is applied. 
    
    EEG = RELAX_filtbutter(EEG, RELAX_cfg.LineNoiseFrequency-3, RELAX_cfg.LineNoiseFrequency+3, 4, 'bandstop' );
    EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, RELAX_cfg.LowPassFilter, 4, 'bandpass' );

    FileName_filt = [FileName_nochan, '-filt'];
    EEG = pop_saveset(EEG, 'filename', FileName_filt, 'filepath', RELAX_cfg.OutputPath);

    %% Apply PREP pipeline functions to detect noisy electrodes.
    %  This deals mainly with flat channels and channels with "improbable"
    %  data.

    noisyOut = findNoisyChannels(EEG);
    
    % Record the electrodes detected as noisy by the PREP findNoisyChannel
    % function.
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
    allchans = {EEG.chanlocs.labels};
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = allchans(1,noisyOut.noisyChannels.all);
    EEG=pop_select(EEG,'nochannel',noisyOut.noisyChannels.all); % Delete noisy electrodes detected by PREP

    %% Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
    %  Epoch data also to detect time periods of noisy data to exclude from MWF cleaning and to delete before wICA cleaning.

    continuousEEG=EEG;
    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(continuousEEG, RELAX_cfg);
    [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg);

    %% Use the continuous data to detect eye-blinks. 
    %  These will be marked in the EEG.event structure as well as in the mask.
    %  The output is continuous data that includes all previous extreme-period markings from the epoched data. 
    %  Remember to have the correct electrodes defined for blink-detection in the RELAX_cfg structure.
    %  Need to exclude the exogenous electrodes from this as they do not have any channel location information.

    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg); % Use an IQR threshold method to detect and mark blinks
        if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0                                       % If a participants doesn't show any blinks, make a note
            NoBlinksDetected{FileNumber,1}=FileName; 
            warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
        end
        if RELAX_cfg.computerawmetrics==1
            [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % Record blink amplitude ratio from raw data for comparison.
        end
    end

end




end