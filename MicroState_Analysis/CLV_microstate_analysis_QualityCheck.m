function CLV_microstate_analysis_QualityCheck()


%% Load in all the data for a give session

session = 'pretest';
reststate = 'restingstate2';

datapath = fullfile(filesep,'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data', session, 'derivatives');
alldirs = {dir(datapath).name};
alldirs_bis = alldirs([dir(datapath).isdir]);
SubDirs = alldirs_bis(contains(alldirs_bis, 'sub'));  % All the subject folders for the current session.

%% To run the Data quality check function, there needs to be a GUI. So the easiest solution is to load the data into an EEGLAB session.

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

for isuj = 1:numel(SubDirs)
    currDir = fullfile(datapath, SubDirs{1,isuj},'restingstate_eeg',filesep);
    currsetAll = {dir(currDir).name};
    currsetIdx = contains(currsetAll, reststate);
    currDataset = currsetAll(currsetIdx); 
    if numel(currDataset)>1
        cs1 = contains(currDataset, '-ssinterp');
        currDataset = currDataset(cs1);
    end
    if ~isempty(currDataset)
        EEG = pop_loadset(currDataset{1,1}, currDir);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
        eeglab redraw
    end

end

setsSelected = 1:numel(ALLEEG);
clusterNumber = 7;

setTables = pop_CheckData(ALLEEG, setsSelected, clusterNumber);






end
