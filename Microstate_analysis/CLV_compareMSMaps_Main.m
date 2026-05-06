function SharedVarTable = CLV_compareMSMaps_Main()

parentDir = '/Users/bolger/Documents/MATLAB/Projects/CeLaVie_EEG';
saveSharedVar = fullfile(parentDir, 'Data', 'Microstate','GrandMean_MSMaps');

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

SharedVar_name = [DataIn_compare.setname,'_shared variances_2.csv'];
SharedVarTable = [];
SharedVarTable = CLV_CompareMSMaps(DataIn_compare, [], 1, 'MetaMap_v2', 6, fullfile(saveSharedVar, SharedVar_name));

end