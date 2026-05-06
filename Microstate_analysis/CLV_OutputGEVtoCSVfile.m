%% Function to extract the GEV from the microstate data of each individual dataset and create a *.csv file that arrranges the data in R long format.
%  Order of the output columns is as follows :
%  Subject number; GEV; Pre-/Post-; Control/Test


% Specify the datasets to load in. 
pretestIn_fpath = fullfile(filesep,'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data','Microstate','MS_SingleSubject_Level','MS_SingleSubject_Sorted',...
    'Sorted_7cluster_solution','Pretest',filesep);
posttestIn_fpath = fullfile(filesep,'Users','bolger','Matlab','Projects','CeLaVie_EEG','Data','Microstate','MS_SingleSubject_Level','MS_SingleSubject_Sorted',...
    'Sorted_7cluster_solution','Posttest',filesep);
% Pretest
pretestdir = dir(pretestIn_fpath);
allpretest_name = {pretestdir.name}; 
sortedpretest = allpretest_name(contains(allpretest_name, 'sortedPublished.set'));  % Cell array of the datasets with sorted microstates.
% Posttest
posttestdir = dir(posttestIn_fpath);
allposttest_name = {posttestdir.name};
sortedposttest = allposttest_name(contains(allposttest_name, 'sortedPublished.set'));

DataIn_pretest = cellfun(@(x) pop_loadset(x,pretestIn_fpath),sortedpretest, 'UniformOutput',false);
DataIn_posttest = cellfun(@(x1) pop_loadset(x1, posttestIn_fpath), sortedposttest, 'UniformOutput',false);

%% Need to check which datasets belong to the Control and Test groups for both pre- and post-test conditions
%  Load in the CeLaVie_SubjectsGroups.xlsx file.

infofileIn = 'CeLaVie_SubjectsGroups.xlsx';
infofileIn_path = '/Users/bolger/Matlab/Projects/CeLaVie_EEG/Data_Processing/';
pretestInfo = readtable([infofileIn_path,infofileIn], ReadVariableNames=true, Sheet='Pretest');
posttestInfo = readtable([infofileIn_path,infofileIn], ReadVariableNames=true, Sheet='Posttest');
% Pretest
codeString_pretest = num2str(pretestInfo.Code);
group_pretest = pretestInfo.Group;
% Posttest
codeString_posttest = cellstr(num2str(posttestInfo.Code));
for cnt = 1:numel(codeString_posttest)
    sptest = isspace(codeString_posttest{cnt,:});
    if sum(sptest)> 0
        codeString_posttest{cnt,1}(find(sptest)) = [];
    end
end
group_posttest = posttestInfo.Group;

currIndx_pretest = zeros(numel(DataIn_pretest),1);
currgroup_pretest = cell(numel(DataIn_pretest),1);
for counter = 1:numel(DataIn_pretest)
    currIndx_pretest(counter,1) = find(cell2mat(cellfun(@(y) contains(DataIn_pretest{1,counter}.setname, y), cellstr(codeString_pretest), 'UniformOutput',false)));
    currgroup_pretest{counter,1} = group_pretest{currIndx_pretest(counter,1),1};
    DataIn_pretest{1, counter}.group = currgroup_pretest{counter,1};
    DataIn_pretest{1, counter}.condition = 'pretest';
end

currIndx_posttest = zeros(numel(DataIn_posttest),1);
currgroup_posttest = cell(numel(DataIn_posttest),1);
for counter2 = 1:numel(currIndx_posttest)
    currIndx_posttest(counter2,1) = find(cell2mat(cellfun(@(y1) contains(DataIn_posttest{1,counter2}.setname, y1),...
        codeString_posttest, 'UniformOutput', false)));
    currgroup_posttest{counter2,1} = group_posttest{currIndx_posttest(counter2,1),1};
    DataIn_posttest{1, counter2}.group = currgroup_posttest{counter2, 1};
    DataIn_posttest{1, counter2}.condition = 'posttest';
end

%% Generate the Control list. 
%  Find the Controls and Test datasets in the Pretest list.

allgroups_pretest = cellfun(@(x) x.group, DataIn_pretest, 'UniformOutput',false);
control_pretest = DataIn_pretest(1,contains(allgroups_pretest,'Control'));   % All control datasets
test_pretest = DataIn_pretest(1,contains(allgroups_pretest, 'Test'));        % ll test datasets.

allgroups_posttest = cellfun(@(x1) x1.group, DataIn_posttest, 'UniformOutput',false);
control_posttest = DataIn_posttest(1,contains(allgroups_posttest, 'Control'));
test_posttest = DataIn_posttest(1,contains(allgroups_posttest,'Test'));



MS4Cluster_nom = {'GEV-4Cluster1';'GEV-4Cluster2';'GEV-4Cluster3';'GEV-4Cluster4'};
MS5Cluster_nom = {'GEV-5Cluster1';'GEV-5Cluster2';'GEV-5Cluster3';'GEV-5Cluster4';'GEV-5Cluster5'}; 
MS6Cluster_nom = {'GEV-6Cluster1';'GEV-6Cluster2';'GEV-6Cluster3';'GEV-6Cluster4';'GEV-6Cluster5';'GEV-6Cluster6'};
MS7Cluster_nom = {'GEV-7Cluster1';'GEV-7Cluster2';'GEV-7Cluster3';'GEV-7Cluster4';'GEV-7Cluster5';'GEV-7Cluster6';'GEV-7Cluster7'};


MS4nom = repmat('MS4', numel(MS4Cluster_nom), 1);
MS5nom = repmat('MS5', numel(MS5Cluster_nom), 1);
MS6nom = repmat('MS6', numel(MS6Cluster_nom), 1);
MS7nom = repmat('MS7', numel(MS7Cluster_nom), 1);

allGEV_noms = cat(1,MS4Cluster_nom, MS5Cluster_nom, MS6Cluster_nom, MS7Cluster_nom);
allMS_noms = cat(1, MS4nom, MS5nom, MS6nom, MS7nom);

%%  Generate tables for Pretest data Controls et Test
% Control-Pretest 
CPretest_new = [];
allctrl_group = repmat('Control', numel(allGEV_noms),1);
allctrl_pretest = repmat('Pretest', numel(allGEV_noms),1);
  
for counter1 = 1:numel(control_pretest)
    S = [];
    ctrlPretest_nom = [];
    ctrl_pretest_MS4Cluster = []; 
    ctrl_pretest_MS5Cluster = [];
    ctrl_pretest_MS6Cluster = [];
    ctrl_pretest_MS7Cluster = [];
    allctrlPretest_control_GEV = [];
    allctrlPretest_nom = [];

    S = split(control_pretest{1,counter1}.setname, '_');
    ctrlPretest_nom = S{1,1};
    ctrl_pretest_MS4Cluster = control_pretest{1,counter1}.msinfo.MSMaps(4).ExpVar';
    ctrl_pretest_MS5Cluster = control_pretest{1,counter1}.msinfo.MSMaps(5).ExpVar';
    ctrl_pretest_MS6Cluster = control_pretest{1,counter1}.msinfo.MSMaps(6).ExpVar';
    ctrl_pretest_MS7Cluster = control_pretest{1,counter1}.msinfo.MSMaps(7).ExpVar';
    allctrlPretest_control_GEV = cat(1, ctrl_pretest_MS4Cluster,ctrl_pretest_MS5Cluster,ctrl_pretest_MS6Cluster, ctrl_pretest_MS7Cluster);  
    allctrlPretest_nom = repmat(ctrlPretest_nom,numel(allGEV_noms),1);

    if contains(control_pretest{1,counter1}.setname, '_restingstate1_')
        restingstate_curr = repmat('Restingstate1', numel(allGEV_noms),1);
    elseif contains(control_pretest{1,counter1}.setname, '_restingstate2_')
        restingstate_curr = repmat('Restingstate2', numel(allGEV_noms),1);
    end
   
    CPretest_curr = cat(2, string(allctrl_group), string(allctrlPretest_nom), string(allctrl_pretest), string(restingstate_curr), allMS_noms, allGEV_noms, allctrlPretest_control_GEV);
    CPretest_new = cat(1,CPretest_new, CPretest_curr);

end

% Test-Pretest 
TPretest_new = [];
alltest_group = repmat('Test', numel(allGEV_noms),1);
alltest_pretest = repmat('Pretest', numel(allGEV_noms),1);
  
for counter1 = 1:numel(test_pretest)
    S = [];
    testPretest_nom = [];
    test_pretest_MS4Cluster = []; 
    test_pretest_MS5Cluster = [];
    test_pretest_MS6Cluster = [];
    test_pretest_MS7Cluster = [];
    alltestPretest_test_GEV = [];
    alltestPretest_nom = [];

    S = split(test_pretest{1,counter1}.setname, '_');
    testPretest_nom = S{1,1};
    test_pretest_MS4Cluster = test_pretest{1,counter1}.msinfo.MSMaps(4).ExpVar';
    test_pretest_MS5Cluster = test_pretest{1,counter1}.msinfo.MSMaps(5).ExpVar';
    test_pretest_MS6Cluster = test_pretest{1,counter1}.msinfo.MSMaps(6).ExpVar';
    test_pretest_MS7Cluster = test_pretest{1,counter1}.msinfo.MSMaps(7).ExpVar';
    alltestPretest_test_GEV = cat(1, test_pretest_MS4Cluster,test_pretest_MS5Cluster,test_pretest_MS6Cluster, test_pretest_MS7Cluster);  
    alltestPretest_nom = repmat(testPretest_nom,numel(allGEV_noms),1);

    if contains(test_pretest{1,counter1}.setname, '_restingstate1_')
        restingstate_curr = repmat('Restingstate1', numel(allGEV_noms),1);
    elseif contains(test_pretest{1,counter1}.setname, '_restingstate2_')
        restingstate_curr = repmat('Restingstate2', numel(allGEV_noms),1);
    end
   
    TPretest_curr = cat(2, string(alltest_group), string(alltestPretest_nom), string(alltest_pretest), string(restingstate_curr), allMS_noms, allGEV_noms, alltestPretest_test_GEV);
    TPretest_new = cat(1,TPretest_new, TPretest_curr);

end

%% need to generate the same as above for the posttest.
%  then save all to a *csv file to load into R.

CPosttest_new = [];
allctrl_group = repmat('Control', numel(allGEV_noms),1);
allctrl_posttest = repmat('Posttest', numel(allGEV_noms),1);

for counter1 = 1:numel(control_posttest)
    S = [];
    ctrlPosttest_nom = [];
    ctrl_posttest_MS4Cluster = []; 
    ctrl_posttest_MS5Cluster = [];
    ctrl_posttest_MS6Cluster = [];
    ctrl_posttest_MS7Cluster = [];
    allctrlPosttest_control_GEV = [];
    allctrlPosttest_nom = [];

    S = split(control_posttest{1,counter1}.setname, '_');
    ctrlPosttest_nom = S{1,1};
    ctrl_posttest_MS4Cluster = control_posttest{1,counter1}.msinfo.MSMaps(4).ExpVar';
    ctrl_posttest_MS5Cluster = control_posttest{1,counter1}.msinfo.MSMaps(5).ExpVar';
    ctrl_posttest_MS6Cluster = control_posttest{1,counter1}.msinfo.MSMaps(6).ExpVar';
    ctrl_posttest_MS7Cluster = control_posttest{1,counter1}.msinfo.MSMaps(7).ExpVar';
    allctrlPosttest_control_GEV = cat(1, ctrl_posttest_MS4Cluster,ctrl_posttest_MS5Cluster,ctrl_posttest_MS6Cluster, ctrl_posttest_MS7Cluster);  
    allctrlPosttest_nom = repmat(ctrlPosttest_nom,numel(allGEV_noms),1);

    if contains(control_posttest{1,counter1}.setname, '_restingstate1_')
        restingstate_curr = repmat('Restingstate1', numel(allGEV_noms),1);
    elseif contains(control_posttest{1,counter1}.setname, '_restingstate2_')
        restingstate_curr = repmat('Restingstate2', numel(allGEV_noms),1);
    end
   
    CPosttest_curr = cat(2, string(allctrl_group), string(allctrlPosttest_nom), string(allctrl_posttest), string(restingstate_curr), allMS_noms, allGEV_noms, allctrlPosttest_control_GEV);
    CPosttest_new = cat(1,CPosttest_new, CPosttest_curr);

end

% Test-Posttest 
TPosttest_new = [];
alltest_group = repmat('Test', numel(allGEV_noms),1);
alltest_posttest = repmat('Posttest', numel(allGEV_noms),1);
  
for counter1 = 1:numel(test_posttest)
    S = [];
    testPosttest_nom = [];
    test_posttest_MS4Cluster = []; 
    test_posttest_MS5Cluster = [];
    test_posttest_MS6Cluster = [];
    test_posttest_MS7Cluster = [];
    alltestPosttest_test_GEV = [];
    alltestPosttest_nom = [];

    S = split(test_posttest{1,counter1}.setname, '_');
    testPosttest_nom = S{1,1};
    test_posttest_MS4Cluster = test_posttest{1,counter1}.msinfo.MSMaps(4).ExpVar';
    test_posttest_MS5Cluster = test_posttest{1,counter1}.msinfo.MSMaps(5).ExpVar';
    test_posttest_MS6Cluster = test_posttest{1,counter1}.msinfo.MSMaps(6).ExpVar';
    test_posttest_MS7Cluster = test_posttest{1,counter1}.msinfo.MSMaps(7).ExpVar';
    alltestPosttest_test_GEV = cat(1, test_posttest_MS4Cluster,test_posttest_MS5Cluster,test_posttest_MS6Cluster, test_posttest_MS7Cluster);  
    alltestPosttest_nom = repmat(testPosttest_nom,numel(allGEV_noms),1);

    if contains(test_posttest{1,counter1}.setname, '_restingstate1_')
        restingstate_curr = repmat('Restingstate1', numel(allGEV_noms),1);
    elseif contains(test_posttest{1,counter1}.setname, '_restingstate2_')
        restingstate_curr = repmat('Restingstate2', numel(allGEV_noms),1);
    end
   
    TPosttest_curr = cat(2, string(alltest_group), string(alltestPosttest_nom), string(alltest_posttest), string(restingstate_curr), allMS_noms, allGEV_noms, alltestPosttest_test_GEV);
    TPosttest_new = cat(1,TPosttest_new, TPosttest_curr);

end

%% Concatenate all elements : Test, Control, Pre- and Post-test lists together and save as a *.csv file. 


MSData_Rlong = 'MSAllData_RLong.csv';
MSData_Rlong_path = '/Users/bolger/Matlab/Projects/CeLaVie_EEG/Data/Microstate';

string2export = cellstr(cat(1, CPretest_new, CPosttest_new, TPretest_new, TPosttest_new));
table2export = cell2table(string2export, 'VariableNames',{'Group', 'Subject', 'Condition', 'Restingstate', 'MSNumber', 'MSClusterIndex', 'GEV'});
writetable(table2export, fullfile(MSData_Rlong_path,MSData_Rlong))







   









