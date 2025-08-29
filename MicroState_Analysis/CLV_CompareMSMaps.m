function sharedVarTable = CLV_CompareMSMaps(EEGIn, IndividualSets, MeanSets, PubSet, Classes, FileName)

%% Verify the input sets.

HasMS = arrayfun(@(x) clv_hasMicrostates(EEGIn(x)), 1:numel(EEGIn));
HasChildren = arrayfun(@(x) clv_DoesItHaveChildren(EEGIn(x)), 1:numel(EEGIn));
HasDyn = arrayfun(@(x) clv_isDynamicsSet(EEGIn(x)), 1:numel(EEGIn));
if ~isempty(PubSet)
    isPublished = arrayfun(@(x) clv_isPublishedSet(EEGIn(x), {PubSet}), 1:numel(EEGIn));
else
    isPublished = zeros(1, numel(EEGIn));
end

AvailableIndSets = find(HasMS & ~HasChildren & ~HasDyn & ~isPublished);     % Finding the individual sets
AvailableMeanSets = find(HasMS & HasChildren & ~HasDyn & ~isPublished);     % Finding the mean sets
AvailablePublishedSets = sum(isPublished);                                  % Finding the published template sets
AvailableSets = [AvailableIndSets, AvailableMeanSets, AvailablePublishedSets];


if isempty(AvailableSets)
    errorMessage = ['No valid sets found for comparing maps. Use ' ...
        '"Tools->Identify microstate maps per dataset" to find and store microstate map data.'];
    if all(matches({'IndividualSets','MeanSets','PublishedSets'}, p.UsingDefaults))
        errorDialog(errorMessage, 'Compare microstate maps error');
        return;
    else
        error(errorMessage);
    end
end

% Validate individual sets
if ~isempty(IndividualSets)
    IndividualSets = unique(IndividualSets, 'stable');
    isValid = ismember(IndividualSets, AvailableIndSets);
    if any(~isValid)
        invalidSetsTxt = sprintf('%i, ', IndividualSets(~isValid));
        invalidSetsTxt = invalidSetsTxt(1:end-2);
        error(['The following individual sets are invalid: ' invalidSetsTxt ...
            '. Make sure you have not selected empty sets, mean sets, dynamics sets, ' ...
            'or sets without microstate maps.']);
    end
end

% Validate mean sets
MeanSetnames = {EEGIn(AvailableMeanSets).setname};
if ~isempty(MeanSets)
    MeanSets = unique(MeanSets, 'stable');
    if isnumeric(MeanSets)
        isValid = ismember(MeanSets, AvailableMeanSets);
        if any(~isValid)
            invalidSetsTxt = sprintf('%i, ', MeanSets(~isValid));
            invalidSetsTxt = invalidSetsTxt(1:end-2);
            error(['The following mean sets are invalid: ' invalidSetsTxt ...
                '. Make sure you have not selected individual datasets, dynamics sets, ' ...
                'or sets without microstate maps.']);
        end
    else
        isValid = ismember(MeanSets, MeanSetnames);
        if any(~isValid)
            invalidSetsTxt = sprintf('%s, ', string(MeanSets(~isValid)));
            invalidSetsTxt = invalidSetsTxt(1:end-2);
            error(['The following mean sets could not be found: ' invalidSetsTxt]);
        else
            % if MeanSets is a string array/cell array of char vectors,
            % convert to integers
            MeanSets = AvailableMeanSets(ismember(MeanSetnames, MeanSets));
        end
    end
end

% Validate PublishedSetsn--- need to redo this one.
if ~isempty(PubSet)
    templatePath = '/Users/bolger/Documents/Documents - MACPRO-BOLGER-01/MATLAB/eeglab2024.0/plugins/MICROSTATELAB2.1/Templates';
    templateNames_all = {dir(templatePath).name};
    if sum(ismember(templateNames_all, [PubSet,'.set'])) == 1
        fprintf('The specified published template set, %s, is in the Templates folder. \n', PubSet)
        PublishedSetIndices = find((ismember(templateNames_all, [PubSet,'.set'])));
    elseif sum(ismember(templateNames_all, [PubSet,'.set'])) == 0
        error(['The following published template sets could not be found in the Templates folder: ' invalidSetsTxt '.\n']);
        PublishedSetIndices = [];
    end

else
    PublishedSetIndices = [];
end


%% Determine whether the comparison is within or across datasets.

% Create a question dialogue box.
SelectedSets = [IndividualSets(:); MeanSets(:)];

if (numel(SelectedSets) + numel(PublishedSetIndices)) == 1
    compWithin = 1;
else
    compWithin = 0;    % Implies comparison across sets
end

SelectedEEG = [];
if ~isempty(SelectedSets)
    SelectedEEG = eeg_store(SelectedEEG, EEGIn(SelectedSets), 1:numel(SelectedSets));
end

nonpublishedSets = 1:numel(SelectedEEG);

if ~isempty(PublishedSetIndices)  % Load in the published template dataset.

    MSTemplate_In = pop_loadset(templateNames_all{PublishedSetIndices}, templatePath);

    SelectedEEG = eeg_store(SelectedEEG, MSTemplate_In, (numel(SelectedEEG)+1):(numel(SelectedEEG)+numel(MSTemplate_In)));
end

%% Check for overlap in cluster solutions

AllMinClasses = arrayfun(@(x) SelectedEEG(x).msinfo.ClustPar.MinClasses, 1:numel(SelectedEEG));
AllMaxClasses = arrayfun(@(x) SelectedEEG(x).msinfo.ClustPar.MaxClasses, 1:numel(SelectedEEG));
MinClasses = max(AllMinClasses);
MaxClasses = min(AllMaxClasses);

if MaxClasses < MinClasses
    errorMessage = 'No overlap in cluster solutions found between all selected sets.';
    if showDlg
        errordlg2(errorMessage, 'Compare microstate maps error');
        return;
    else
        error(errorMessage);
    end
end

classes = MinClasses:MaxClasses;
classChoices = sprintf('%i Classes|', classes);

%% Check for different channel numbers and create common set of channels
if numel(SelectedEEG) > 1
    nChannels = arrayfun(@(x) numel(SelectedEEG(x).chanlocs), 1:numel(SelectedEEG));
    [~, minSetIdx] = min(nChannels);

    for i=1:numel(SelectedEEG)
        if i == minSetIdx
            continue
        end

        [LocalToGlobal, ~] = MakeResampleMatrices(SelectedEEG(i).chanlocs, SelectedEEG(minSetIdx).chanlocs);

        for class=SelectedEEG(i).msinfo.ClustPar.MinClasses:SelectedEEG(i).msinfo.ClustPar.MaxClasses
            SelectedEEG(i).msinfo.MSMaps(class).Maps = SelectedEEG(i).msinfo.MSMaps(class).Maps*LocalToGlobal';
            SelectedEEG(i).chanlocs = SelectedEEG(minSetIdx).chanlocs;
        end
    end
end

%% Calculate and export shared variances.

if compWithin
    nClasses = 0;
end

MapCollection    = [];
CLabelCollection = [];
if compWithin
    for i = MinClasses:MaxClasses
        MapCollection = [MapCollection; SelectedEEG.msinfo.MSMaps(i).Maps];
        for j = 1:i
            CLabelCollection  = [CLabelCollection,sprintf("%s (%i)",SelectedEEG.msinfo.MSMaps(i).Labels{j},i)];
        end
    end
else
    for i=1:numel(SelectedEEG)
        MapCollection = [MapCollection; SelectedEEG(i).msinfo.MSMaps(Classes).Maps];
        for j=1:Classes
            CLabelCollection = [CLabelCollection, sprintf("%s (%s)", SelectedEEG(i).msinfo.MSMaps(Classes).Labels{j}, SelectedEEG(i).setname)];
        end
    end
end
CorrMat = MyCorr(double(MapCollection)').^2;
sharedVarTable = array2table(CorrMat * 100,'VariableNames',CLabelCollection,'RowNames',CLabelCollection);  % Write the covariance matrix to a table. 
writetable(sharedVarTable, FileName, 'WriteRowNames',true);

assignin('base', 'SelectedEEG', SelectedEEG);
assignin('base', 'Classes', Classes);
assignin('base', 'FileName', FileName);

FileName = Compare_DisplayMDS_MSSolutions(SelectedEEG, Classes, FileName);  % Call of function to display the grandmean MS maps versus template maps on a multidimensional scale. 



end