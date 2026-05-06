% EEGLAB history file generated on the 21-Nov-2025
% ------------------------------------------------

fnameAll = {};
fpathAll = {};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

cd("/Users/bolger/Documents/MATLAB/Projects/CeLaVie_EEG/");
data_fpath = fullfile(cd, 'Data','Posttest', 'Derivatives');
ls = dir(data_fpath);
allfolders = {ls.name};
issubject = contains(allfolders, "sub");
currfolders = allfolders(issubject);


for subcount = 1:numel(currfolders)

    dircurr = fullfile(data_fpath, currfolders{1,subcount}, 'restingstate_eeg');
    isset = contains({dir(dircurr).name},'ssinterp.set');
    dircurr_names = {dir(dircurr).name};
    currfolders2load = dircurr_names(1,isset);
    allcurrdir = repmat({dircurr},1,size(currfolders2load,2));
    DataIn = cellfun(@(x,y) pop_loadset(x,y), currfolders2load, allcurrdir, 'UniformOutput',false); % Load in the datasets.
    ALLEEG = cat(2, ALLEEG, DataIn{:});

end


[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
eeglab redraw;

%% Calculate individual-level microstate maps for each condition seperately (Pre-test and Post-test)

ClustPar.UseAAHC = 0;
ClustPar.MinClasses = 4;
ClustPar.MaxClasses = 7;
ClustPar.Restarts = 20;
ClustPar.MaxMaps = inf;
ClustPar.GFPPeaks = 1;
ClustPar.IgnorePolarity = 1;
ClustPar.Normalize = 1;

[ALLEEG, CurrentSet, com] = pop_FindMSMaps(ALLEEG, 1:numel(ALLEEG), "ClustPar", ClustPar);

savepath = '/Users/bolger/Documents/MATLAB/Projects/CeLaVie_EEG/Data/Microstate/MS_Posttest_SingleSubject_Level/Posttest_microstates';
savefiles = {EEG(:).setname};

for scnt = 1:numel(ALLEEG)
    [EEGout, ~] = pop_saveset(ALLEEG(scnt), [ALLEEG(scnt).setname,'_microstates.set'], savepath);
    eeglab redraw;
end

%% Combine the individual-level microstates to create a grandmean microstata set.

testfiles = {'MS_Pretest_SingleSubject_Level','MS_Posttest_SingleSubject_Level'};
dataFpath_ms = fullfile(cd, 'Data','Microstate');
allfpaths_ms = cellfun(@(x2) fullfile(dataFpath_ms, x2), testfiles, 'UniformOutput',false);
allfiles_ms = cellfun(@(x3) dir(x3), allfpaths_ms, 'UniformOutput',false);
AllEEGIn = [];

for scnt = 1:numel(allfiles_ms)

    allfilesms_names = {allfiles_ms{1,scnt}.name};
    isset = find(contains(allfilesms_names, '.set'));
    AllFiles_curr = char(allfilesms_names(isset));
    alleeg_curr = pop_loadset('filename', cellstr(AllFiles_curr), 'filepath', allfpaths_ms{1,scnt});
    AllEEGIn = cat(1, AllEEGIn, alleeg_curr');

end

SelectedSets = 1:numel(AllEEGIn);
meanName = 'Grandmean_microstate';
EEG_gm = pop_CombMSMaps(AllEEGIn,SelectedSets , 'MeanName', meanName,'IgnorePolarity', true);

savepath_gm = '/Users/bolger/Matlab/Projects/CeLaVie_EEG/Data/Microstate/GrandMean_MSMaps';
[EEGout_gm, ~] = pop_saveset(EEG_gm, meanName, savepath_gm);

%% Call of function CLV_sortMSMaps_Main() sort the individual level microstate clusters based on the sorted grand-mean microstate cluster.

[EEGSorted, Path2Figure] = CLV_sortMSMaps_Main();

%% Call of function to compare the grandmean microstate map against the template microstate map.
%  This function can also be used to compare inidividual microstate
%  maps against the grandmean template map.

SharedVarTable = CLV_compareMSMaps_main();

%% Call of function to look for outliers at the individual subject level. 

fnameAll2test = {};
fpathAll2test = {};
fnameAll = {};
fpathAll = {};

[fname2test, fpath2test] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets to sort ', 'multiselect', 'on');
[parentdir, ~,~] = fileparts(fpath2test);
[~,~,~] = fileparts(fname2test);
main_dir = fileparts(parentdir);
fnameAll2test = [fnameAll2test; {fname2test}];
fpathAll2test = [fpathAll2test; {fpath2test}];

if ~ischar(fname2test)                                   % If only loading in one dataset.
    fprintf('Loading in the datasets to be tested for outliers...\n');
    fpathAll2test = repmat(fpathAll2test,1,size(fname2test,2));   % Add repeated copies of filepath to have same dimensions as filenames.
    DataIn_test = cellfun(@(x,y) pop_loadset(x,y), fnameAll2test{1,1}, fpathAll2test, 'UniformOutput',false); % Load in the datasets.
    DataIn_test = cell2mat(DataIn_test);
elseif ischar(fname2compare)
    fprintf('Loading in the dataset to be tested for outliers...\n')
    DataIn_test = pop_loadset(fnameAll2test, fpathAll2test{1,1});
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
pop_DetectOutliers(DataIn_test, 1:numel(DataIn_test), 'Classes', 6)

%% Call of function to backfit the microstate maps to the EEG.
%  Here we use the grandmean template to backfit.
%  Load in the data to be re-expressed and the grandmean template to be
%  applied for the backfitting.

fnameAllFitMS = {};
fpathAllFitMS = {};

[fname2fit, fpath2fit] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets to be fitted ', 'multiselect', 'on');
[fnametemplate, fpathtemplate] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the template dataset for backfitting ', 'multiselect', 'on');
fnameAllFitMS = [fnameAllFitMS; {fname2fit}];
fpathAllFitMS = [fpathAllFitMS; {fpath2fit}];
fpathAllFitMS = repmat(fpathAllFitMS, 1, numel(fname2fit));  % Generate a file path for each dataset to be fit.

fprintf('Loading in the datasets to be backfitted...\n');
DataIn_FitMS = cellfun(@(x,y) pop_loadset(x,y), fnameAllFitMS{1,1}, fpathAllFitMS, 'UniformOutput',false);
fprintf('Laoding in the template dataset to apply for backfitting...\n');
DataIn_FitTemp = pop_loadset(fnametemplate, fpathtemplate);
DataIn_FitMS = cell2mat(DataIn_FitMS);
fields2add_temp = setdiff(fieldnames(DataIn_FitMS(1)), fieldnames(DataIn_FitTemp)); % Find the fields that are missing from the DataIn_FitTemp structure.
for cnt = 1:numel(fields2add_temp)
    DataIn_FitTemp.(fields2add_temp{cnt,1}) = [];
end

if numel(fieldnames(DataIn_FitTemp))==numel(fieldnames(DataIn_FitMS))
    fprintf('The inidividual microstate datasets and the template dataset have the same number of fields.\n');
    DataIn_BackFitAll = [DataIn_FitMS, DataIn_FitTemp];
else
    fprintf('The inidividual microstate datasets and the template dataset do not have the same number of fields.\n')
end

% Create structure to define the backfitting parameters.
FitParams = [];
FitParams.PeakFit = 1; % Back fit based on GFP peaks.
FitParams.b = 0;       % This is ignored if fitting based on the GFP peaks.
FitParams.lambda = 0; % The penality for non-smoothness in temporal smooting but ignored if fitting based on the GFP peaks.
% Note that we can backfit based on each individual dataset's own
% microstates or based on the grandmean template. For the moment try both.

[EEG, ~, ~] = pop_FitMSMaps(DataIn_BackFitAll, 1:numel(DataIn_BackFitAll)-1, 'FitPar', FitParams, 'TemplateSet', numel(DataIn_BackFitAll));

savepath = '/Users/bolger/Documents/MATLAB/Projects/CeLaVie_EEG/Data/Microstate/MS_Posttest_SingleSubject_Level/Posttest_backfitted_GrandMean_metatemp_6Class';
EEG = arrayfun(@(s1) pop_saveset(s1, s1.filename, savepath, 'savemode','onefile'), EEG, 'UniformOutput', false);

%% Export the backfitted temporal parameters for analysis in Ragu.
fnameAllSaveMS = {};
fpathAllSaveMS = {};

[fname2save, fpath2save] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets whose microstate parameters you wish to save for further analysis ', 'multiselect', 'on');
fnameAllSaveMS = [fnameAllSaveMS; {fname2save}];
fpathAllSaveMS = [fpathAllSaveMS; {fpath2save}];
fpathAllSaveMS = repmat(fpathAllSaveMS, 1, numel(fname2save));  % Generate a file path for each dataset to be saved.
% Load in actual EEG datasets.
fprintf('Loading in the datasets to be saved for further analysis...\n');
DataIn_SaveMS = cellfun(@(x,y) pop_loadset(x,y), fnameAllSaveMS{1,1}, fpathAllSaveMS, 'UniformOutput',false);

setsSelected = 1:numel(DataIn_SaveMS);
s = split(fpath2save, '/');
saveFileName = [s{end-1},'_MSParameters_6cluster_backfitGrandMean_export.csv'];
SaveFilePath = fullfile(fpath2save,saveFileName);

MS_stats = pop_SaveMSParameters(cell2mat(DataIn_SaveMS), setsSelected, 'Classes', 6, 'Filename', SaveFilePath);

%% 


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

fnameIn1 = {};
fpathIn1 = {};
[fname2load, fpath2load] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets ', 'multiselect', 'on');
fnameIn1 = [fnameIn1; {fname2load}];
fpathIn1 = [fpathIn1; {fpath2load}];
fpathIn1 = repmat(fpathIn1, 1, numel(fname2load));
fprintf('Loading in the first set of datasets for analysis in Ragu...\n')
dataIn_1 = cellfun(@(x,y) pop_loadset(x,y), fnameIn1{1,1}, fpathIn1, 'UniformOutput',false);
dataIn_1 = cell2mat(dataIn_1);

fnameIn2 = {};
fpathIn2 = {};
[fname2load1, fpath2load1] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets ', 'multiselect', 'on');
fnameIn2 = [fnameIn2; {fname2load1}];
fpathIn2 = [fpathIn2; {fpath2load1}];
fpathIn2 = repmat(fpathIn2, 1, numel(fname2load1));
fprintf('Loading in the second set of datasets for analysis in Ragu...\n')
dataIn_2 = cellfun(@(x,y) pop_loadset(x,y), fnameIn2{1,1}, fpathIn2, 'UniformOutput',false);
dataIn_2 = cell2mat(dataIn_2);

% Concatenate the datasets from both conditions.
dataInAll = [dataIn_1, dataIn_2];

RScurr = '_restingstate2_';
allnoms = {dataInAll(:).setname};
RScurr_indx = find(contains(allnoms, RScurr));

ALLEEG = dataInAll(RScurr_indx);
eeglab redraw



