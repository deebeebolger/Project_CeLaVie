%% Function to store data in BIDS structure.
%  This script takes a directory of raw files from either the pre-test or
%  the post-test and creates a directory in BIDS structure. The raw files
%  are renamed according to BIDS convention and moved to the new BIDS
%  structured directory.
%
%  Overview of the structure of the filenames:
%  sub-XXX_ses-pretest_task-restingstate_eeg.bdf'
%  There needs to be a folder structure like:
%  ses-pretest -- sub-XXX -- restingstate1 -- eeg
%  ses-pretest -- sub-XXX -- restingstate2 -- eeg
%  ses-posttest -- sub-XXX -- restingstate1 -- eeg
%  ses-posttest -- sub-XXX -- restingstate2 -- eeg
%  File structure:
%  data.file : name of file with data path
%  data.restingstate  ==> 1 or two
%  data.session     ==> pretest = 1 posttest = 2
%

%% Set up raw data files : each file will contain restingstate1 and restingstate2
%  ----------------------------------------------------------------------------

datapath = fullfile(filesep, 'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data',filesep);   % Create BIDS file structure here.
raw_datapath = fullfile(filesep, 'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data_RestingState','RELAXRaw', filesep);

files2load = {dir(fullfile(raw_datapath,'*.bdf')).name}; % Return names of files with *.bdf extension as cell array.
session = 'posttest';
datatype = 'eeg';

%% Create the BIDS based on each raw file.

for fcount = 1:numel(files2load)
    
    currfile = files2load{1,fcount};
    currfile_split = strsplit(currfile, '_');
    x = strsplit(currfile_split{1,2}, '.');
    sujnum_curr = x{1,1};
    subject = ['sub-0',sujnum_curr];

    if strcmp(currfile_split{1,1}, 'RS1') 
        currtask = 'restingstate1';
    elseif strcmp(currfile_split, 'RS2')
        currtask = 'restingstate2';
    end
    datapath_bids = fullfile(datapath, subject, currtask, datatype);
    mkdir(datapath_bids)

    newname = [subject,'_ses-',session,'_task-',currtask,'_',datatype,'.',x{1,2}];
    movefile(fullfile(raw_datapath,currfile), fullfile(datapath_bids, newname));

end

