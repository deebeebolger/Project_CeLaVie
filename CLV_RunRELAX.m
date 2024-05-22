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

    %% RESAMPLE THE DATA HERE. 

    newFS = 512;
    fprintf('Downsampling from %d to %d Hz', EEG.srate, newFS);
    EEG = pop_resample(EEG, newFS);

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
            NoBlinksDetected{fcounter,1}=FileName;
            warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
        end
        if RELAX_cfg.computerawmetrics==1
            [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % Record blink amplitude ratio from raw data for comparison.
        end
    end

    %% Record extreme artifact rejection details for all participants in single table.

    RELAXProcessingExtremeRejectionsAllParticipants(fcounter,:) = struct2table(epochedEEG.RELAXProcessingExtremeRejections,'AsArray',true);
    rawEEG=continuousEEG;          % Take a copy of the not yet cleaned data for calculation of all cleaning SER and ARR at the end

    %% Mark artifacts for calculating SER and ARR, regardless of whether MWF is performed (RELAX v1.1.3 update).

    if RELAX_cfg.computecleanedmetrics==1 && (RELAX_cfg.Do_MWF_Once==0 || RELAX_cfg.Do_MWF_Twice==0 || RELAX_cfg.Do_MWF_Thrice==0)

        [Marking_artifacts_for_SER_ARR, ~] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1;
        [Marking_artifacts_for_SER_ARR] = RELAX_horizontaleye(continuousEEG, RELAX_cfg);
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1;
        [Marking_artifacts_for_SER_ARR, ~] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1;
        Marking_all_artifacts_for_SER_ARR.RELAX.NaNsForExtremeOutlierPeriods=continuousEEG.RELAX.NaNsForExtremeOutlierPeriods;
        [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out.
        
        if isfield(continuousEEG.RELAX,'eyeblinkmask')
            Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(continuousEEG.RELAX.eyeblinkmask==1)=1;
            [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'blinks');
        end
        
        continuousEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength;
        rawEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength;
    else
        fprintf('No need to run this marking of artifacts for SER and ARR calculation as: \n Do MWF Once is %d\n Do MWF Twice is %d\n Do MWF Trice is %d.\n',...
            RELAX_cfg.Do_MWF_Once, RELAX_cfg.Do_MWF_Twice, RELAX_cfg.Do_MWF_Thrice)

    end

    %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK AND CLEAN MUSCLE ARTIFACTS. Carry out first round of MWF. 
    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask  
    
     if RELAX_cfg.Do_MWF_Once==1

        %% Use epoched data and FFT to detect slope of log frequency log power, add periods exceeding muscle threshold to mask.
        
        [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);  
        if RELAX_cfg.computerawmetrics==1
            [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg); % record muscle contamination metrics from raw data for comparison.
        end

        EEG=continuousEEG; % Return continuousEEG to the "EEG" variable for MWF processing

        %% If including eye blink cleaning in first round MWF, then insert eye blink mask into noise mask.
        
        if RELAX_cfg.MWFRoundToCleanBlinks==1

            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))= NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');

        else
            fprintf('Not including eye blink cleaning in this first round.\nThis will be done in round %d.\n', RELAX_cfg.MWFRoundToCleanBlinks)
        end

        %% The following pads very brief lengths of mask periods in the template (without doing this, very short periods can lead to rank deficiency), 
        % and excludes extreme artifacts from the cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the clean or artifact mask, but does apply cleaning to them).
        
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out.
        
        EEG.RELAX.NoiseMaskFullLengthR1=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR1=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        
        [EEG] = RELAX_perform_MWF_cleaning(EEG, RELAX_cfg);          

        EEG.RELAXProcessingRoundOne=EEG.RELAXProcessing; % Record MWF cleaning details from round 1 in EEG file          
        RELAXProcessingRoundOne=EEG.RELAXProcessingRoundOne; % Record MWF cleaning details from round 1 into file for all participants
        
        if isfield(RELAXProcessingRoundOne,'Details')
            RELAXProcessingRoundOne=rmfield(RELAXProcessingRoundOne,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundOne,'Details')
                EEG.RELAXProcessingRoundOne=rmfield(EEG.RELAXProcessingRoundOne,'Details');
            end
        end
        
        %% Record processing statistics for all participants in single table:
        
        RELAXProcessingRoundOneAllParticipants(fcounter,:) = struct2table(RELAXProcessingRoundOne,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        
        % Save round 1 MWF pre-processing:
        if RELAX_cfg.saveround1==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF'])
            end
            SaveSetMWF1 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF', filesep FileName '_MWF1.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF1 ); 
        end

     end  % End of Do MWF Once if condition.
     
     %% PERFORM A SECOND ROUND OF MWF. THIS IS HELPFUL IF THE FIRST ROUND DOESN'T SUFFICIENTLY CLEAN ARTIFACTS. 

     % This has been suggested to be useful by Somers et al (2018)
     % (particularly when used in a cascading fashion). 

     % However, I can see risks. If artifact masks fall on task relevant
     % activity in both rounds of the MWF, it may be that the task relevant data
     % is just cleaned right out of the signal.
     
     if RELAX_cfg.Do_MWF_Twice==1

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        
        %% If blinks weren't initially detected because they were  disguised by the the muscle artifact, detect them here (this happens in <1/200 cases, but is a good back up).
        
        if RELAX_cfg.ProbabilityDataHasNoBlinks==0
            if EEG.RELAX.IQRmethodDetectedBlinks(1,1)==0
                continuousEEG=EEG;
                [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
                EEG=continuousEEG;
            end
        end

        %% If including eye blink cleaning in second round MWF, then insert eye blink mask into noise mask:
        
        if isfield(EEG.RELAX, 'eyeblinkmask')
            if RELAX_cfg.MWFRoundToCleanBlinks==2
                fprintf('Cleaning eye-blinks in this second round of MWF.\n');
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
            end
        end

        %% The following pads very brief lengths of mask periods in the template (without doing this, very short periods can lead to rank deficiency),
        %  and excludes extreme artifacts from the cleaning template (so the MWF cleaning step just ignores extreme artifacts in it's template - 
        %  doesn't include them in either the clean or artifact mask, but does apply cleaning to them).
        
        [EEG] = RELAX_pad_brief_mask_periods(EEG, RELAX_cfg, 'blinks');
        
        EEG.RELAX.NoiseMaskFullLengthR2=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        
        [EEG] = RELAX_perform_MWF_cleaning(EEG, RELAX_cfg);           
 
        EEG.RELAXProcessingRoundTwo=EEG.RELAXProcessing;     % Record MWF cleaning details from round 2 in EEG file
        RELAXProcessingRoundTwo=EEG.RELAXProcessingRoundTwo; % Record MWF cleaning details from round 2 into file for all participants
        
        if isfield(RELAXProcessingRoundTwo,'Details')
            RELAXProcessingRoundTwo=rmfield(RELAXProcessingRoundTwo,'Details');
        end

        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundTwo,'Details')
                EEG.RELAXProcessingRoundTwo=rmfield(EEG.RELAXProcessingRoundTwo,'Details');
            end
        end

        %% Record processing statistics for all participants in single table:
        
        RELAXProcessingRoundTwoAllParticipants(fcounter,:) = struct2table(RELAXProcessingRoundTwo,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        
        % Save round 2 MWF pre-processing:
        if RELAX_cfg.saveround2==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF'])
            end
            SaveSetMWF2 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF', filesep FileName '_MWF2.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF2 ); 
        end     

     end

     %% PERFORM A THIRD ROUND OF MWF. 

     if RELAX_cfg.Do_MWF_Thrice == 1

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;

        %% If isfield(EEG.RELAX,'ProportionMarkedInMWFArtifactMaskTotalR2') 
        %  NWB added to make sure function doesn't bug when trying to check this variable if it doesn't exist

            if EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2<0.05
                if isfield(EEG.RELAX, 'eyeblinkmask')
                    EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                    EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                    EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
                end
            end

        %% Epoch the data into 1 second epochs with a 500ms overlap. 
        %  Outputs both the ContinuousEEG (which has been filtered above by this point) and the epoched data as EEG.
        
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);

        %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK ARTIFACTS

        [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask
        
        % Use the filtered continuous data to detect horizontal eye movements and mark these in the EEG.event as well as in the mask.
        % You may want to simply reject horizontal eye movements at a later stage if your task requires participants to look straight ahead
        % for the entire task. 
        % Alternatively, if your task requires participants to complete horizontal eye movements time locked to
        % a stimuli, this section will mark every event with these horizontal eye movements as an artifact, and should not be implemented.
        
        % The output is continuous data:
        [continuousEEG] = RELAX_horizontaleye(continuousEEG, RELAX_cfg);

        %% Return to the "EEG" variable for MWF processing:
        
        EEG=continuousEEG;
        
        % If including eye blink cleaning in third round MWF, then insert eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==3
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
        end

        %% The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks');
        
        EEG.RELAX.NoiseMaskFullLengthR3=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR3=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);               
        
        if isfield(EEG.RELAX, 'eyeblinkmask')            % if eyeblinkmask has been created, do the following.
            EEG.RELAX=rmfield(EEG.RELAX,'eyeblinkmask'); % remove variables that are no longer necessary
        end
        
        EEG.RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 in EEG file
        RELAXProcessingRoundThree=EEG.RELAXProcessing;     % Record MWF cleaning details from round 3 into file for all participants
        
        if isfield(RELAXProcessingRoundThree,'Details')
            RELAXProcessingRoundThree=rmfield(RELAXProcessingRoundThree,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundThree,'Details')
                EEG.RELAXProcessingRoundThree=rmfield(EEG.RELAXProcessingRoundThree,'Details');
            end
        end

        %% Record processing statistics for all participants in single table:
        
        RELAXProcessingRoundThreeAllParticipants(fcounter,:) = struct2table(RELAXProcessingRoundThree,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');

        if RELAX_cfg.saveround3==1
            if ~exist([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '3xMWF'], 'dir')
                mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '3xMWF'])
            end
            SaveSetMWF3 =[RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep '3xMWF', filesep FileName '_MWF3.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF3 ); 
        end

     end % End of Do MWF Thrice if condition.

     %% Perform robust average re-referencing of the data and reject periods marked as extreme outliers.

     if RELAX_cfg.Do_MWF_Once==0
        EEG=continuousEEG;
     end
     [EEG] = RELAX_average_rereference(EEG);
     EEG = eeg_checkset( EEG );


    %% Reject periods that were marked as NaNs in the MWF masks because they showed extreme shift within the epoch or extremely improbable data.

    EEG = eeg_eegrej(EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion); % Takes everything out for subject 18!!

    %% Perform wICA on ICLabel identified artifacts that remain.
    %  The following performs wICA, implemented on only the components marked as artifact by ICLabel.
    % fastica_symm is repeated up to 3 times in the case of non-convergence to ensure non-convergence doesn't impair cleaning. 
    % fastica performs quickly on continuous data and doesn't seem to be inferior at cleaning compared to AMICA (which is much slower). 
    % It also seems to be comparable (or only slightly worse) than extended infomax (run via cudaICA for speed).
    
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
       
        EEG.RELAXProcessing_wICA.aFileName=cellstr(FileName);

        [EEG,~, ~, ~, ~] = RELAX_wICA_on_ICLabel_artifacts(EEG,RELAX_cfg.ICA_method, 1, 0, EEG.srate, 5,'coif5',RELAX_cfg.Report_all_ICA_info); 
        
        % adding 'Report_all_wICA_info' to the end of the parameters specified will optionally report proportion of ICs categorized as each category, 
        % and variance explained by ICs from each category (function is ~20s slower if this is implemented)
        EEG = eeg_checkset( EEG );

        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(fcounter,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end

    %% Perform ICA subtract on ICLabel identified artifacts that remain:
    %  The following performs ICA sutraction, implemented on only the components marked as artifact by ICLabel.
    
    if RELAX_cfg.Perform_ICA_subtract==1
    
        EEG.RELAXProcessing_ICA.aFileName=cellstr(FileName);
        EEG = RELAX_ICA_subtract(EEG,RELAX_cfg);
        EEG = eeg_checkset( EEG );

        RELAXProcessing_ICA=EEG.RELAXProcessing_ICA;
        
        % Record processing statistics for all participants in single table:
        RELAXProcessing_ICA_AllParticipants(fcounter,:) = struct2table(RELAXProcessing_ICA,'AsArray',true);
    end
    
    EEG.RELAX.Data_has_been_cleaned=1;

    %% COMPUTE CLEANED METRICS.

    if RELAX_cfg.computecleanedmetrics==1    
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        [continuousEEG, ~] = RELAX_metrics_blinks(continuousEEG, epochedEEG);
        [continuousEEG, ~] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg);

        [continuousEEG] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG); % this is only a good metric for testing only the cleaning of artifacts marked for cleaning by MWF, see notes in function.

        EEG=continuousEEG;
        EEG = rmfield(EEG,'RELAXProcessing');  % More weird housekeeping...

        if isfield(EEG,'RELAX_Metrics')
            if isfield(EEG.RELAX_Metrics, 'Cleaned')
                if isfield(EEG.RELAX_Metrics.Cleaned,'BlinkAmplitudeRatio')
                    CleanedMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
                    CleanedMetrics.BlinkAmplitudeRatio(CleanedMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Cleaned,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
                if isfield(EEG.RELAX_Metrics.Cleaned,'All_SER')
                    CleanedMetrics.All_SER(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_SER;
                    CleanedMetrics.All_ARR(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_ARR;
                end
            end
            if isfield(EEG.RELAX_Metrics, 'Raw')
                if isfield(EEG.RELAX_Metrics.Raw,'BlinkAmplitudeRatio')
                    RawMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio;
                    RawMetrics.BlinkAmplitudeRatio(RawMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Raw,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    RawMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    RawMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
            end   
        end
    end

    %% Record warnings about potential issues:
    
    EEG.RELAX_issues_to_check.aFileName=cellstr(FileName);
    if size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)>RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1); % 1.1.4: fix dimension specification error
    else
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=0;
    end
    if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1))...
            >=RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
            (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1));
    else
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=0;
    end
    if EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier>0.20
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier;
    else 
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=0;
    end
    if isfield(EEG.RELAX, 'IQRmethodDetectedBlinks') % if IQRmethodDetectedBlinks has been created, do the following (thanks to Jane Tan for the suggested bug fix when IQRmethodDetectedBlinks is not created)
        EEG.RELAX_issues_to_check.NoBlinksDetected=(EEG.RELAX.IQRmethodDetectedBlinks==0);
    end
    if RELAX_cfg.Do_MWF_Once==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R1=isa(EEG.RELAXProcessingRoundOne.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Twice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R2=isa(EEG.RELAXProcessingRoundTwo.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R3=isa(EEG.RELAXProcessingRoundThree.RankDeficiency,'char');
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        if EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge(1,3);
    end
    if RELAX_cfg.Perform_ICA_subtract==1
        if EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_ICA.fastica_symm_Didnt_Converge(1,3);
    end
    
    %% If defined, interpolate the rejected electrodes after cleaning.

    if strcmp(RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning,'yes')
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    end



 end




end