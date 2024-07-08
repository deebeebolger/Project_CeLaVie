%% Function to create the RELAX configuration structure.

function Relax_cfg = CLV_CreateRELAX(testtype)

    Relax_cfg = [];  % Initialise the RELAX configuration structure.
    Relax_cfg.caploc = fullfile(filesep,'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data_Processing','ChanLocs128.mat');  % Path to change. Path to channel location file.
    Relax_cfg.myPath = fullfile(filesep,'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data'); % Path to change. Path to the data to be processed.
    Relax_cfg.myPathRaw = fullfile(Relax_cfg.myPath, testtype);                           % Path to change. Path to raw data folder. All datasets here will be processed.
    Relax_cfg.OutputPath = fullfile(Relax_cfg.myPath, 'RELAXProcessed');                  % Path to change. Path to folder in which processed datasets are saved. This may need to be created.
    
    %% Set the initial extreme outlier cleaning parameters.
    Relax_cfg.sample_rate = 1024;
    Relax_cfg.ms_per_sample = (1/Relax_cfg.sample_rate)*1000;

    Relax_cfg.Do_MWF_Once = 1;
    Relax_cfg.Do_MWF_Twice = 1;
    Relax_cfg.Do_MWF_Thrice = 1;
    
    Relax_cfg.Perform_wICA_on_ICLabel = 1; 
    Relax_cfg.Perform_ICA_subtract = 1;
    Relax_cfg.ICA_method = 'infomax';
    Relax_cfg.Report_all_ICA_info = 'no';
    
    Relax_cfg.computerawmetrics = 1;                            % Compute blink and muscle metrics from the raw data.
    Relax_cfg.computecleanedmetrics = 1;                        % Compute SER, ARR, blink and muscle metrics from the cleaned data.
    Relax_cfg.MWFRoundToCleanBlinks = 2;                        % Clean blinks on first, second round. 
    Relax_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks = 'no'; % Otherwise too restrictive
    
    Relax_cfg.ProbabilityDataHasNoBlinks = 0;                   % 0 = data almost certainly has blinks, 1 = data might not have blinks, 2 = data definitely doesn't have blinks.
    Relax_cfg.DriftSeverityThreshold = 10;                      % MAD from the median of all electrodes. This could be set lower and would catch less severe drift.
    Relax_cfg.ProportionWorstEpochsForDrift = 0.30;             % Maximum proportion of epochs to include in the mask from drift artifact type.
    
    Relax_cfg.ExtremeVoltageShiftThreshold = 20;                % Threshold MAD from the median all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe voltage shifts within the epoch.
    Relax_cfg.ExtremeAbsoluteVoltageThreshold = 500;            % Microvolts max or min above which will be excluded from cleaning and deleted from data
    Relax_cfg.ExtremeImprobableVoltageDistributionThreshold = 8;% Threshold SD from the mean of all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe improbable data.
    Relax_cfg.ExtremeSingleChannelKurtosisThreshold = 8;        % Threshold kurtosis of each electrode against the same electrode in different epochs. This could be set lower and would catch less severe kurtosis.
    Relax_cfg.ExtremeAllChannelKurtosisThreshold = 8;           % Threshold kurtosis across all electrodes. This could be set lower and would catch less severe kurtosis.
    Relax_cfg.ExtremeDriftSlopeThreshold = -4;                  % Slope of log frequency log power below which to reject as drift without neural activity.
    Relax_cfg.ExtremeBlinkShiftThreshold = 8;                   % How many MAD from the median across blink affected epochs to exclude as extreme data.

    Relax_cfg.MinimumArtifactDuration = 1200;                   % In ms. It's better to make this value longer than 1000ms, as doing so will catch diminishing artifacts that aren't detected in a neighbouring 1000ms period, which might still be bad.
    Relax_cfg.MinimumBlinkArtifactDuration = 800;               % Blink marking is based on the maximum point of the blink rather than the 1000ms divisions for muscle artifacts, so this can be shorter than the value above (blinks do not typically last >500ms).
    Relax_cfg.BlinkElectrodes = {'C29';'C17';'C16';'C30';'C8';'C28';'C27';'C21';'C13';'C10'};
    Relax_cfg.HEOGLeftpattern = ["C30", "D7", "D8", "D9", "D23", "D10", "D22", "D24", "C31"]; % Lateral electrodes maximally effected by hEOGs.
    Relax_cfg.HEOGRightpattern = ["C8", "C7","B27","B28","B26", "B29", "B24", "B14", "C9"];   % Right lateral electrodes maximally affected by hEOGs.
    Relax_cfg.BlinkMaskFocus = 150;                             % This value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window.  
    Relax_cfg.HorizontalEyeMovementType = 2;                    % 1 to use the IQR method, 2 to use the MAD method for identifying threshold. IQR method less effective for smaller sample sizes (shorter files).
    Relax_cfg.HorizontalEyeMovementThreshold = 2;               % MAD deviation from the median that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
    Relax_cfg.HorizontalEyeMovementThresholdIQR = 1.5;          % If IQR method set above, IQR deviation that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
    Relax_cfg.HorizontalEyeMovementTimepointsExceedingThreshold = 25; % The number of timepoints (ms) that exceed the horizontal eye movement threshold within the test period (set below) before the period is marked as horizontal eye movement.
    Relax_cfg.HorizontalEyeMovementTimepointsTestWindow = (2*Relax_cfg.HorizontalEyeMovementTimepointsExceedingThreshold)-1;
    Relax_cfg.HorizontalEyeMovementFocus = 200;                 % Buffer window, masking periods earlier and later than the time where horizontal eye movements exceed the threshold.
    
    Relax_cfg.HighPassFilter = 0.25;                            % Sets the high pass filter. 1Hz is best for ICA decomposition if you're examining just oscillatory data, 0.25Hz seems to be the highest before ERPs are adversely affected by filtering.
    Relax_cfg.LowPassFilter = 80;
    Relax_cfg.LineNoiseFrequency = 50;
    Relax_cfg.ElectrodesToDelete = {'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp', 'EXG1','EXG2', 'EXG3',...
        'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8'};
    Relax_cfg.KeepAllInfo = 0;                                  % Setting this value to 1 keeps all the details from the MWF pre-processing and MWF computation. Helpful for debugging if necessary but makes for large file sizes.
    Relax_cfg.saveextremesrejected = 0;                         % Setting this value to 1 tells the script to save the data after only filtering, extreme channels have been rejected and extreme periods have been noted.
    
    Relax_cfg.saveround1 = 1;
    Relax_cfg.saveround2 = 1;
    Relax_cfg.saveround3 = 1;
    Relax_cfg.OnlyIncludeTaskRelatedEpochs = 0; 

    Relax_cfg.MuscleSlopeThreshold = -0.59;                     % Log-frequency log-power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle.
    Relax_cfg.MaxProportionOfDataCanBeMarkedAsMuscle = 0.50;    % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
    Relax_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 0.05;
    Relax_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.05;
    Relax_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.20;
    Relax_cfg.InterpolateRejectedElectrodesAfterCleaning = 'yes';
    Relax_cfg.MWFDelayPeriod_ms = 40;                         % Define the MWF delay period in ms.
    Relax_cfg.MWFDelayPeriod = floor((Relax_cfg.MWFDelayPeriod_ms*(1/1000))*Relax_cfg.sample_rate); % The MWF includes both spatial and temporal information when filtering out artifacts. Longer delays apparently improve performance. 
    
  
end 