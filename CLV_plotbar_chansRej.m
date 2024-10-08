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

%% Load in the datasets selected above and extract the number of rejected channels from all datasets.

dataIn = cellfun(@(x,y) pop_loadset(x,y), fnameAll, fpathAll, 'UniformOutput',false); 
rejChans = cellfun(@(x1) x1.RELAXProcessingExtremeRejections.PREPBasedChannelToReject, dataIn, 'UniformOutput',false);
subjectTags = cellfun(@(x2) x2.setname(1:7), dataIn, 'UniformOutput',false);
chanrejNum = cellfun(@(x3) numel(x3), rejChans, 'UniformOutput',false);
SRate = cellfun(@(s1) s1.srate, dataIn, 'UniformOutput', false);

%% Extract the total duration of rejected time for all datasets.
%  Then calculate the total duration for each dataset. 

rejTIntvals = cellfun(@(y1) y1.RELAX.ExtremelyBadPeriodsForDeletion, dataIn, 'UniformOutput', false);
rejTDur = cellfun(@(y2, s2) diff(y2')/s2, rejTIntvals, SRate, 'UniformOutput',false);
rejTDur_total = cellfun(@(y3) sum(y3), rejTDur, 'UniformOutput',false);

%% Plot bar-chart of the rejected channels for the selected participants.

chanrejBar_fig = figure;
set(chanrejBar_fig, 'Color', [1 1 1], 'Name', 'Bar chart of rejected channels per subject');

X = categorical(subjectTags);
Y = cell2mat(chanrejNum);
b = bar(X,Y);
b.BarWidth = 0.45;
b.FaceColor = [0 .5 .5];
b.EdgeColor = [0 .9 .9];
b.LineWidth = 1.5;
labelsBar = string(b.YData);
xtips = b.XEndPoints;
ytips = b.YEndPoints;
text(xtips, ytips, labelsBar, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

%% Plot bar-chart of the rejected time intervals for each of the selected participants.

timerejBar_fig = figure;
set(timerejBar_fig, 'Color', [1 1 1], 'Name', 'Bar chart of total rejected time per subject.')

X1 = categorical(subjectTags);
Y1 = cell2mat(rejTDur_total);
b1 = bar(X1, Y1);
b1.BarWidth = 0.45;
b1.FaceColor = [0 .5 .5];
b1.EdgeColor = [0 .9 .9];
b1.LineWidth = 1.5;
labelsBar1 = string(b1.YData);
xtips1 = b1.XEndPoints;
ytips1 = b1.YEndPoints;
text(xtips1, ytips1, labelsBar1, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

end 



