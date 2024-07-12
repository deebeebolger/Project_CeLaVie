%% Function to store data in BIDS structure.
%  This script takes a directory of raw files from either the pre-test or
%  the post-test and creates a directory in BIDS structure. The raw files
%  are renamed according to BIDS convention and moved to the new BIDS
%  structured directory.
% 
%  Each subject's data corresponds to a directory of raw data containing
%  subdirectories for each session and data modality. 
%  For each subject, there is a metadata file, "dataset_description_eeg.json", 
%  file with details of the experiment task and the EEG recording system
%  and experimental setup.
%
%  
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
tasktype = 'restingstate';

%% 

xlspath = '/Users/bolger/Documents/Projects/CeLaVie/Celavie_docs';
xlsfname = 'Pre-Post-Test_EEG_Celavie_metadata.xlsx';
MetaDataIn = CLV_LoadMetaData(xlspath, xlsfname);

%% Create participants.json file.

Age = MetaDataIn.Age;
Handedness = MetaDataIn.MainDom_;
Gender = MetaDataIn.Genre;
School = MetaDataIn.Ecole;
SujetCode = MetaDataIn.Code;     % Subject Code
Sujetnumber = MetaDataIn.Sujet;  % Subject Number

Handedness_u = unique(Handedness);
mtHandedness = cellfun(@isempty, Handedness_u,'UniformOutput', false);
HandednessU = Handedness_u(~cell2mat(mtHandedness)); 

Gender_u = unique(Gender);
mtGender = cellfun(@isempty, Gender_u, 'UniformOutput', false);
GenderU = Gender_u(~cell2mat(mtGender));

School_u = unique(School);
mtSchool = cellfun(@isempty, School_u, 'UniformOutput', false);
SchoolU = School_u(~cell2mat(mtSchool));

participants = struct();  % Initialise empty structure;
participants.subjectnum.Description = 'the number of the participant';
participants.subjectnum.Type = 'ordinal';

participants.subjectcode.Description = 'the code attributed to each participant';
participants.subjectcode.Type = 'anonymous code';

participants.age.Description = "age of the participant";
participants.age.Units = "years";

participants.handedness.Description = 'handedness of the participant as reported by the participant.';
[participants.handedness.Levels.right participants.handedness.Levels.both participants.handedness.Levels.left] = deal(HandednessU{:});

participants.gender.Description = 'gender of the participant';
[participants.gender.Levels.female participants.gender.Levels.male] = deal(GenderU{:});

participants.school.Description = 'school of the participant';
[participants.school.Levels.AL participants.school.Levels.M] = deal(SchoolU{:});

% Save general participants information structure as *.json
json_participants = jsonencode(participants, PrettyPrint=true);
json_participants_title = 'participants.json';
fid = fopen(fullfile(datapath,session, json_participants_title), 'w');
fprintf(fid, '%s', json_participants);
fclose(fid);

%% Create the file participants.tsv
%  Need to extract the following information from the MetaData excel file:
%  - subject number
%  - age
%  - handedness
%  - gender
%  - school

% Extract the anonymisation codes from individual dataset titles.
txt = files2load;
pat = digitsPattern;
subCodes = extract(txt, pat);
startzero_ind = find(startsWith(subCodes, '0'));
CodesOnly = cellfun(@(x) x(2:end), subCodes(startzero_ind), 'UniformOutput',false);

sujetsIndex = find(ismember(SujetCode, str2double(CodesOnly)));
participants_info = T(sujetsIndex, ["Code","Sujet","Age","Genre","MainDom_", "Ecole"]);
participants_csv_path = fullfile(datapath, session, 'participants.csv');
writetable(participants_info,participants_csv_path)

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
    datapath_bids = fullfile(datapath, session, subject, tasktype, datatype, currtask);
    mkdir(datapath_bids)

    newname = [subject,'_ses-',session,'_task-',tasktype,'_',datatype,'_',currtask.',x{1,2}];
    movefile(fullfile(raw_datapath,currfile), fullfile(datapath_bids, newname));

end

% Create derivatives directory
derivpath_bids = fullfile(datapath, session, 'derivatives');

