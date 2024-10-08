
function [filesinBIDS2Load, pathinBIDS,  derivpathinBIDS, allsubjects_title, session_path, participants_info] = CLV_createBIDS(session, datatype, tasktype)
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

datasetInfo_path = fullfile(filesep, 'Users','bolger','Matlab','Projects','CeLaVie_EEG', 'Data_Processing');
datasetInfo_fname = 'subjectlevel_eeg.csv';


%% CHECK IF A SESSION RELATED FOLDER EXISTS AND IF NOT, CREATE.

session_path = fullfile(datapath, session);
if exist("session_path","dir")
    fprintf('The session path : %s already exists.\n', session_path)
elseif ~exist("session_path","dir")
    sess_status = mkdir(session_path);     % Generate the path.
end

%% CREATE A DERIVATIVES DIRECTORY IN WHICH TO SAVE PROCESSED DATA.

derivpath_bids = fullfile(session_path, 'derivatives');

%% LOAD IN THE META DATA FROM EXCEL FILE.

xlspath = '/Users/bolger/Documents/Projects/CeLaVie/Celavie_docs';
xlsfname = 'Pre-Post-Test_EEG_Celavie_metadata.xlsx';
MetaDataIn = CLV_LoadMetaData(xlspath, xlsfname);
assignin('base', "MetaDataIn", MetaDataIn)

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
fid = fopen(fullfile(session_path, json_participants_title), 'w');
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

sujetsIndex = find(ismember(SujetCode, str2double(subCodes)));
participants_info = MetaDataIn(sujetsIndex, ["Code","Sujet","Age","Genre","MainDom_", "Ecole", "Diametre", "AxeAntero_Posterieur",...
    "AxeGauche_Droite", "Couleur", "Taille"]);
participants_csv_path = fullfile(session_path, 'participants.csv');
if ~exist(participants_csv_path, 'file')
    fprintf('A participants file: \n %s already exists, \n create a new participant info file. \n', participants_csv_path);
    writetable(participants_info,participants_csv_path)
else
    fprintf('A participants file: \n %s already exists,\n add current participant info to existing file. \n',participants_csv_path);
    T = readtable(participants_csv_path);
    TAll = [T; participants_info];
    writetable(TAll, participants_csv_path)
end

%% Create the BIDS based on each raw file.
filesinBIDS2Load = cell(numel(files2load),1);
pathinBIDS = cell(numel(files2load),1);
derivpathinBIDS = cell(numel(files2load),1);
allsubjects_title = cell(numel(files2load),1);

for fcount = 1:numel(files2load)

    currfile = files2load{1,fcount};
    currfile_split1 = strsplit(currfile, '_');
    x = strsplit(currfile_split1{1,2}, '.');
    pat1 = digitsPattern;
    xdigit = extract(currfile_split1{1,1}, pat1);
    if length(xdigit{1,1})==2
        subject = ['sub-0',xdigit{1,1}];
    elseif length(xdigit{1,1})==1
        subject = ['sub-00',xdigit{1,1}];
    end
    allsubjects_title{fcount,1} = subject;

    if strcmp(x{1,1}, 'RS1')
        currtask = 'restingstate1';
    elseif strcmp(x{1,1}, 'RS2')
        currtask = 'restingstate2';
    end
    datapath_bids = fullfile(session_path, subject, tasktype, datatype, currtask);  % Contains raw data.
    if ~exist(datapath_bids, "dir")
        status1 = mkdir(datapath_bids);
    elseif exist(datapath_bids, "dir")
        fprintf('The bids datapath, %s, already exists.\n', datapath_bids);
    end

    derivpath_curr = fullfile(derivpath_bids, subject, [tasktype,'_',datatype]);    % Contains processed data and related files.
    if ~exist(derivpath_curr, "dir")
        status2 = mkdir(derivpath_curr);
    elseif exist(derivpath_curr, "dir")
        fprintf('The bids datapath, %s, already exists.\n', derivpath_curr);
    end

    newname = [subject,'_sess-',session,'_task-',tasktype,'_',datatype,'_',currtask,'.',x{1,2}];
    movefile(fullfile(raw_datapath,currfile), fullfile(datapath_bids, newname));
    
    filesinBIDS2Load{fcount, 1} = newname;
    currsujet_info = MetaDataIn(sujetsIndex(fcount), ["Diametre", "AxeAntero_Posterieur","AxeGauche_Droite", "CommentairesElectrodesRS1"]);
    make_eeg_json(datasetInfo_path, datasetInfo_fname, newname(1:end-4), datapath_bids, currsujet_info);

    assignin('base', 'currsujet_info', currsujet_info)

    pathinBIDS{fcount, 1} = datapath_bids;
    derivpathinBIDS{fcount, 1} = derivpath_curr;

end
end % end of function

function make_eeg_json(datapath, datafname, currfname, eegdatapath, currsujet_info)

    currpath = fullfile(datapath, datafname);
    T = readtable(currpath, "ReadRowNames",true, "Delimiter",',', 'ReadVariableNames',false);
    
    Tbis = rows2vars(T);
    Tbis.HeadCircumference = [num2str(currsujet_info.Diametre),'cm'];
    Tbis.HeadAntero_Posterior = [num2str(currsujet_info.AxeAntero_Posterieur), 'cm'];
    Tbis.HeadLeft_Right = [num2str(currsujet_info.AxeGauche_Droite), 'cm'];
    Tbis.SubjectArtefactDescription = [currsujet_info.CommentairesElectrodesRS1];
    
    json_eeg = jsonencode(table2struct(Tbis), "PrettyPrint",true);
    json_eeg_title = [currfname,'_eeg.json'];
    fid = fopen(fullfile(eegdatapath, json_eeg_title), 'w');
    fprintf(fid, '%s', json_eeg);
    fclose(fid);

end 



