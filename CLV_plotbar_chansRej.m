function CLV_plotbar_chansRej()
% Function to create a bar chart presenting the rejected channels par
% participant.
% For this routine, the user needs to open an EEGLAB session and load in
% the datasets containing information on the rejected channels (index and
% channel names).
% *************************************************************************


%% Select the datasets whose rejected channels you wish to present.
rechoix = 1;
prompt = {'Load another dataset (Yes/No):'};

fnameAll = {};
fpathAll = {};

while rechoix == 1

    [filename, filepath] = uigetfile({'*.bdf; *.set' 'BDF or eeglab set file';'*.bdf' 'BDF'; '*.set' 'eeglab set file'}, 'Choose a BDF or .set data file ', 'multiselect', 'on');
    [parentdir, ~,~] = fileparts(filepath);
    [~,name,ext] = fileparts(filename);
    main_dir = fileparts(parentdir);
    fnameAll = [fnameAll; {filename}];
    fpathAll = [fpathAll; {filepath}];
    
    dlgtitle = 'Input';
    fieldsize = [1 45];
    definput = {'Yes'};
    answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

    if strcmp(answer, 'Yes')
        rechoix = 1;
    else
        rechoix = 0;
    end
    
end 

%% Load in the datasets selected above.

dataIn = cellfun(@(x,y) pop_loadset(x,y), fnameAll, fpathAll, 'UniformOutput',false); 
rejChans = cellfun(@(x1) x1.RELAXProcessingExtremeRejections.PREPBasedChannelToReject, dataIn, 'UniformOutput',false);
subjectTags = cellfun(@(x2) x2.setname(1:7), dataIn, 'UniformOutput',false);
chanrejNum = cellfun(@(x3) numel(x3), rejChans, 'UniformOutput',false);

%% Plot bar-chart of the rejected channels for the selected participants.

chanrejBar_fig = figure;
set(chanrejBar_fig, 'Color', [1 1 1], 'Name', 'Bar chart of rejected channels per subject');

X = categorical(subjectTags);
Y = cell2mat(chanrejNum);
b = bar(X,Y)
b.BarWidth = 0.45;
b.FaceColor = [0 .5 .5];
b.EdgeColor = [0 .9 .9];
b.LineWidth = 1.5;
labelsBar = string(b.YData);
xtips = b.XEndPoints;
ytips = b.YEndPoints;
text(xtips, ytips, labelsBar, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')





