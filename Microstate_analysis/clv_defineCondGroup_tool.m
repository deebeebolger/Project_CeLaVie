%% Function to define the condition and group of each individual dataset.
function clv_defineCondGroup_tool()

condition = 'Posttest';
fnameIn = {};
fpathIn = {};

[fname2corr, fpath2corr] = uigetfile({ '*.set' 'choose .set files'}, 'Choose the datasets to be corrected for group and condition ', 'multiselect', 'on');
fnameIn = [fnameIn; {fname2corr}];
fpathIn = [fpathIn; {fpath2corr}];
fpathIn = repmat(fpathIn, 1, numel(fname2corr));  % Generate a file path for each dataset to be fit.

fprintf('Loading in the %s datasets to be corrected...\n', condition);
DataIn_corr = cellfun(@(x,y) pop_loadset(x,y), fnameIn{1,1}, fpathIn, 'UniformOutput',false);

% Need to load in the excel file defining the participant groups.
InfoFile_path = '/Users/bolger/Documents/MATLAB/Projects/CeLaVie_EEG/Data_Processing/';
InfoFile_name = 'CeLaVie_SubjectsGroups.xlsx';
InfoFile_fullpath = fullfile(InfoFile_path, InfoFile_name);

DatasetInfo = readtable(InfoFile_fullpath, "FileType","spreadsheet", "Sheet",condition);

for counter = 1:numel(DataIn_corr)
    namecurr = DataIn_corr{1,counter}.setname;
    S = split(namecurr,'-');
    S1 = S{2,1}(1:3);
    if strcmp(S1(1), '0')
        Scurr = S1(2:3);
    else
        Scurr = S1;
    end
    currCode = str2double(Scurr);
    codeIndx = find([DatasetInfo.Code == currCode]);
    currGroup = DatasetInfo.Group(codeIndx);
    currGroup = currGroup{1,1};
    DataIn_corr{1,counter}.group = currGroup;
    DataIn_corr{1,counter}.condition = condition;
    DataIn_corr{1,counter}.subject = Scurr;

    EEGcurr = pop_saveset(DataIn_corr{1,counter}, DataIn_corr{1,counter}.filename, DataIn_corr{1,counter}.filepath, 'savemode', 'onefile');

end

end 