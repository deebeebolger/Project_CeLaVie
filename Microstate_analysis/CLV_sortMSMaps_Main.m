function [EEGsorted, fig_fullPath] = CLV_sortMSMaps_Main()

%% Get sorting parameters.

promptIn = {'Template dataset type (published (publised) or grandmean (GMdataset))', 'Number of cluster classes to sort (all for all clusters)', 'Template no. of classes', 'Ignore polarity'};
dlgtitle = 'Define sorting parameters.';
sizeoffield = [1 40; 1 40; 1 40; 1 40];
defaultInputs = {'published', 'all', '7', '1'};
paramsIn = inputdlg(promptIn, dlgtitle, sizeoffield, defaultInputs);

templateDataset = paramsIn{1,1};
ClassNumber = paramsIn{2,1};
templateClassNumber = str2double(paramsIn{3,1});
PolarityIgnore = str2double(paramsIn{4,1});

%% Define paths to the data to be sorted.

fnameAll2sort = {};
fpathAll2sort = {};
fnameAll = {};
fpathAll = {};

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

%% Specifying the template set type and loading in the template data.

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

end


EEGsorted = CLV_sortMSMaps(DataIn_All, selectedSets, templateSet, ClassNumber, templateClassNumber, PolarityIgnore, 0, 0);

    %% Show the sorted MS maps
    ShowMaps = 1;
    nameSplit = cellfun(@(b) split(b,"_"),{EEGsorted.setname}, 'UniformOutput',false);
    savepath = {EEGsorted.filepath};
    classes = [EEGsorted(1).msinfo.ClustPar.MinClasses:EEGsorted(1).msinfo.ClustPar.MaxClasses];

    if ShowMaps

        fig_h = arrayfun(@(a) pop_ShowIndMSMaps(a, 1:numel(a), 'Classes', classes, 'Visible', 0), EEGsorted, 'UniformOutput', false);
        figname  = cellfun(@(x1) [x1,'_sorted_metatemplate_v2.fig'],{EEGsorted.setname}, 'UniformOutput', false);
        fig_fullPath = cellfun(@(x,y) fullfile(x,y), savepath, figname, 'UniformOutput', false);
        cellfun(@(x,y) saveas(x,y), fig_h, fig_fullPath, 'UniformOutput', false);
    end

end