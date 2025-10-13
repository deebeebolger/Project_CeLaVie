function CLV_microstate_analysis()
%% ***************************************************************************************
% Date : Juin 2025    Programmed by : D. BOLGER
% Script to carry out the microstate analysis of the CeLaVie data using the
% Microstatelab toolbox, a EEGLAB plugin (Kalburgi et al, 2024).
% The script carries out the following steps in the analysis of microstates
% of spontaneous data.
%
%
%% Load in the datasets

fnameAll = {};
fpathAll = {};

[filename, filepath] = uigetfile({'*.bdf; *.set' 'BDF or eeglab set file';'*.bdf' 'BDF'; '*.set' 'eeglab set file'}, 'Choose a BDF or .set data file ', 'multiselect', 'on');
[parentdir, ~,~] = fileparts(filepath);
[~,~,~] = fileparts(filename);
main_dir = fileparts(parentdir);
fnameAll = [fnameAll; {filename}];
fpathAll = [fpathAll; {filepath}];
analysisType = 'compare_Microstates';    % Or 'find_Microstates'  'sortMaps' 'combine_Microstates'


if ~ischar(filename)

    fpathAll = repmat(fpathAll,1,size(filename,2));   % Add repeated copies of filepath to have same dimensions as filenames.
    DataIn = cellfun(@(x,y) pop_loadset(x,y), fnameAll{1,1}, fpathAll, 'UniformOutput',false); % Load in the datasets.
    DataIn = cell2mat(DataIn);
elseif ischar(filename)
    DataIn = pop_loadset(fnameAll, fpathAll{1,1});
end

allsetnames = {DataIn.setname};
setname_split = split(allsetnames, '_'); 
allconditions = unique(setname_split(:,:,1));
if numel(allconditions)>1
    conditionCurr = [];
end


if strcmp(analysisType, 'find_Microstates')

    fprintf('------- Carrying out first-level micro-state analysis-----------\n')


    %% Carry out the microstate segmentation of the individual datasets.

    SelectedSets = 1:numel(DataIn);   % Indices of the datasets on which to carry out MS segmentation.
    ClustPar = [];                  % Cluster paramters to define and call in the pop_FindMSMaps () function.
    ClustPar.UseAAHC = false;       % Carries out kmeans, if true it carries out AAHC
    ClustPar.MinClasses = 4;        % Minimum number of clusters to identify
    ClustPar.MaxClasses = 7;        % Maximum number of clusters to identify
    ClustPar.Restarts = 20;           % Number of times kmeans algorithm is restarted (this is ignored if using AAHC)
    ClustPar.MaxMaps = inf;         % Maximum number of data samples to use to identify clusters.
    ClustPar.GFPPeaks = true;       % Whether clustering should be limited to GFP peaks.
    ClustPar.IgnorePolarity = true; % whether maps with inverted polarity should be considered part of the same cluster.
    ClustPar.Normalize = true;      % set to false if using AAHC.
    TTFrD = false;

    DataOut = [];
    savepath = fullfile(cd, 'Data', 'Microstate');

    DataOut = CLV_findMS_maps(DataIn, SelectedSets, ClustPar, TTFrD);
 
    for counter = 1:numel(SelectedSets)
        currfilename = [DataOut(counter).setname, '_MSMaps'];
        EEG = pop_saveset(DataOut(counter), 'filename', currfilename, 'filepath', savepath);
    end

elseif strcmp(analysisType, 'combine_Microstates')

    %% Add restingstate, condition and group information to all the datasets.
    % Start by detecting controls and tests (group-level)
    conds = {'Pretest', 'Posttest'};

    % Load in excel file with participant-group information.
    % We are interested in the Code and Group columns.
    excelpath = '/Users/bolger/Matlab/Projects/CeLaVie_EEG/Data_Processing';
    xls_fname = 'CeLaVie_SubjectsGroups.xlsx';

    for counter = 1:numel(conds)
        
        TIn = readtable(fullfile(excelpath, xls_fname), 'Sheet', conds{1,counter});  
        condcurr = [lower(conds{1,counter}(1)),conds{1,counter}(2:end)];
        for counter2 = 1:numel(TIn.Code)
            currcode = TIn.Code(counter2);
            Idx = find(contains({DataIn.setname},['-0',num2str(currcode),'_']));
            Idxcond = contains({DataIn(Idx).setname}, condcurr);
            IdxCurr = Idx(Idxcond);
            for icnt = 1:numel(IdxCurr)
                DataIn(IdxCurr(icnt)).group = TIn.Group(counter2);
                DataIn(IdxCurr(icnt)).condition = condcurr;
            end  
        end
    end
    
    %% Now define the group, condition and resting-state blocks across which you want to calculate the mean.
    %  The programme will look to the group and condition fields of each
    %  dataset.

    AverageAcross = {{'Test','Control'}, {'pretest','posttest'}, {'restingstate1_', 'restingstate2'}};    % Averaging group, condition and restingstate criteria.
  
    Addbin = contains([DataIn.group], AverageAcross{1,1}) + contains({DataIn.condition}, AverageAcross{1,2})+contains({DataIn.setname}, AverageAcross{1,3}); % Finding the datasets that include all average criteria.
    selectedSets = find(Addbin==3);
    CurrentDatasets = {DataIn(selectedSets).setname}';

    fprintf('Average across the following datasets : \n')
    cellfun(@(x) fprintf('- %s\n',x), CurrentDatasets)

    cd('/Users/bolger/Matlab/Projects/CeLaVie_EEG')
    MSFig_savepath = fullfile(cd, 'Data', 'Microstate', 'Figures');
    meanNom = ['Mean_',AverageAcross{1,1},'_',AverageAcross{1,2},'_',AverageAcross{1,3}];  % Define a title for the current mean MS dataset.

    if iscell(meanNom)
        meanNom = strcat(meanNom{1,:});
    end

    EEGOUT = [];
    EEGOUT = CLV_combineMSMaps(DataIn, selectedSets, meanNom, 1, 1, MSFig_savepath);       % Call of function to find average over individual MS maps.

    %% Save the EEGOUT structure that contains the mean microstate data; the number of microstates ranges from the minimum to the
    %  the maximum number of microstates (generally 4 to 7).

    savepath = fullfile(cd, 'Data', 'Microstate');
    currfilename = [meanNom, '_MSMaps'];
    EEG = pop_saveset(EEGOUT, 'filename', currfilename, 'filepath', savepath);

elseif strcmp(analysisType, 'sortMaps')

    %% Load in the datasets to be sorted as well as any template maps, other than published template maps.

    fnameAll2sort = {};
    fpathAll2sort = {};
    fnameAll = {};
    fpathAll = {};

    templateDataset = 'GMdataset'; % Could also 'GMdataset' (Grand Mean dataset) or 'own'

    [fname2sort, fpath2sort] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets to sort ', 'multiselect', 'on');
    [parentdir, ~,~] = fileparts(fpath2sort);
    [~,~,~] = fileparts(fname2sort);
    main_dir = fileparts(parentdir);
    fnameAll2sort = [fnameAll2sort; {fname2sort}];
    fpathAll2sort = [fpathAll2sort; {fpath2sort}];

    %% Load in the actual data.
    if ~ischar(fname2sort)
        fprintf('Loading in the datasets to be sorted...\n');
        fpathAll2sort = repmat(fpathAll2sort,1,size(fname2sort,2));   % Add repeated copies of filepath to have same dimensions as filenames.
        DataIn_sort = cellfun(@(x,y) pop_loadset(x,y), fnameAll2sort{1,1}, fpathAll2sort, 'UniformOutput',false); % Load in the datasets.
        DataIn_sort = cell2mat(DataIn_sort);
    elseif ischar(fname2sort)
        fprintf('Loading in the dataset to be sorted...\n')
        DataIn_sort = pop_loadset(fnameAll2sort, fpathAll2sort{1,1});
    end

    %% Specifying the templaetset type and loading in the template data. 

    if strcmp(templateDataset, 'published')
        % open dialogue box to choose the published dataset to select as
        % the template set.
        %set up a uigetfile to pick out the publised dataset that you wish
        %to use as template.
        templateSet = templateDataset;
        [MSMappub_name, MSMappub_path] = uigetfile({'*.set' 'choose a .m file'}, 'Choose a published template dataset', 'multiSelect', 'off');
        fprintf('The following published templateset will be used as the template dataset : %s. \n', MSMappub_name);
        DataTempPub = pop_loadset(MSMappub_name, MSMappub_path);
        DataIn_All = {DataIn_sort, DataTempPub};
        selectedSets = 1:numel(DataIn_All{1,1});
        templateSet = DataTempPub.setname;

    elseif strcmp(templateDataset, 'GMdataset')

        [fname2temp, fpath2temp] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the template dataset ', 'multiselect', 'off');  % Should only select one.
        fprintf('The following Grand-mean dataset will be used as the template set : %s. \n', fname2temp);

        fprintf('Loading in the template dataset...\n');
        DataIn_temp = pop_loadset(fname2temp, fpath2temp);
        DataIn_All = {DataIn_sort, DataIn_temp};
        selectedSets = 1:numel(DataIn_All{1,1});
        templateSet  = 'GMdataset';

    elseif strcmp(templateDataset, 'own')


    end


    ClassNumber = 'all';
    templateClassNumber = 7;
    PolarityIgnore = 1;

    EEGOutput = CLV_sortMSMaps(DataIn_All, selectedSets, templateSet, ClassNumber, templateClassNumber, PolarityIgnore, 0, 0); % Call of function to sort the grandmean microstate maps according to a published template.

    %% Show the sorted MS maps
    ShowMaps = 1;
    nameSplit = cellfun(@(b) split(b,"_"),{EEGOutput.setname}, 'UniformOutput',false);
    savepath = fullfile(cd, 'Data', 'Microstate', 'Figures');
    classes = [EEGOutput(1).msinfo.ClustPar.MinClasses:EEGOutput(1).msinfo.ClustPar.MaxClasses];

    if ShowMaps

        fig_h = arrayfun(@(a) pop_ShowIndMSMaps(a, 1:numel(a), 'Classes', classes, 'Visible', 0), EEGOutput, 'UniformOutput', false);
        savepath = {EEGOutput.filepath};
        figname  = {EEGOutput.setname};
        fig_fullPath = cellfun(@(x,y) fullfile(x,y), savepath, figname, 'UniformOutput', false);
        cellfun(@(x,y) saveas(x,y), fig_h, fig_fullPath, 'UniformOutput', false);
    end

elseif strcmp(analysisType, 'compareMaps')
    
    parentDir = '/Users/bolger/Matlab/Projects/CeLaVie_EEG/';
    saveSharedVar = fullfile(parentDir, 'Data', 'Microstate','GrandMean_MSMaps');

    %% Load in the datasets to be compared.

    fnameAll2compare = {};
    fpathAll2compare = {};
    fnameAll = {};
    fpathAll = {};

    [fname2compare, fpath2compare] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets to sort ', 'multiselect', 'on');
    [parentdir, ~,~] = fileparts(fpath2compare);
    [~,~,~] = fileparts(fname2compare);
    main_dir = fileparts(parentdir);
    fnameAll2compare = [fnameAll2compare; {fname2compare}];
    fpathAll2compare = [fpathAll2compare; {fpath2compare}];

    %% Load in the actual data.

    if ~ischar(fname2compare)                                   % If only loading in one dataset. 
        fprintf('Loading in the datasets to be sorted...\n');
        fpathAll2compare = repmat(fpathAll2compare,1,size(fname2compare,2));   % Add repeated copies of filepath to have same dimensions as filenames.
        DataIn_compare = cellfun(@(x,y) pop_loadset(x,y), fnameAll2compare{1,1}, fpathAll2sort, 'UniformOutput',false); % Load in the datasets.
        DataIn_compare = cell2mat(DataIn_compare);
    elseif ischar(fname2compare)
        fprintf('Loading in the dataset to be compared...\n')
        DataIn_compare = pop_loadset(fnameAll2compare, fpathAll2compare{1,1});
    end

    SharedVar_name = [DataIn_compare.setname,'_shared variances.csv'];
    SharedVarTable = [];
    SharedVarTable = CLV_CompareMSMaps(DataIn_compare, [], 1, 'Custo2017', 4, fullfile(saveSharedVar, SharedVar_name));

end


end


