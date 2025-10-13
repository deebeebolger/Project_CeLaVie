function EEGout = CLV_sortMSMaps(EEGIn, selectSets, TemplateSet, ClassNum, TemplateClassNum, IgnorePolarity, varargin)
% Function to sort the micorstate maps based on either :
% - A published microstate template set
% - An already sorted grandmean microstate map, which is used as a template
% to sort its children.
% It based on the pop_SortMSMaps function form
% MICROSTATELAB: The EEGLAB toolbox for resting-state microstate analysis
% Version 1.0
% Authors:
% Thomas Koenig (thomas.koenig@upd.unibe.ch)
% Delara Aryan  (dearyan@chla.usc.edu)

% The input parameters :
% - EEGIn : all the EEG datasets loaded into the session (structure)
%
% - selectSets : the indices of those datasets on which to carry out
% sorting (vector)
%
% - TemplateSet : specifies the template set to use for sorting. Can either
% be "own", one of the published templates such as Custo2017, "manual" for
% manual sorting or an index to specify one of the datasets of EEGIn.
%
% - ClassNum : the class number/s indicating which cluster solutions to
% sort. The default is "all". If the TemplateSet is "manual" you can only
% specify a single class.
%
% - TemplateClassNum : a scalar indicating the cluster solution to apply
% when sorting. If TemplateSet is "manual" this can be ignored.
%
% - IgnorePolarity : Consider maps with inverted polarities as the same
% class. Default is 1.
%
% Optiobal inputs :
% - Stepwise : Used only when sorting using "own" maps to apply stepwise
% sorting to sort all class solutions.
%
% - SortOrder : This is only applied in the case of manual sorting. Array of microstate map indices
% to apply for manual reordering.
%
% - NewLabel : This is only applied for manual sorting. A string or cell
% array of character vectors defining new microstate map labels.

%% Extract the optional parameter values.

assignin('base', 'varargin', varargin)

if length(varargin)==1
    Stepwise = varargin{1,1};   % Defines whether or not to carry out stepwise sorting. Default is 0.
    if strcmp(TemplateSet, 'manual')
        SortOrder = 1;
        NewLabel  = 1;
    else
        SortOrder = 0;
        NewLabel  = 0;
    end
elseif length(varargin) == 2
    Stepwise = varargin{1,1};
    SortOrder = varargin{1,2};
    if strcmp(TemplateSet, 'manual')
        NewLabel  = 1;
    else
        NewLabel  = 0;
    end
elseif length(varargin) == 3
    Stepwise = varargin{1,1};
    SortOrder = varargin{1,2};
    NewLabel = varargin{1,3};
elseif isempty(varargin)

    if strcmp(TemplateSet, 'manual')
        Stepwise  = 1;
        SortOrder = 1;
        NewLabel  = 1;
    else
        Stepwise  = 0;
        SortOrder = 0;
        NewLabel  = 0;
    end

elseif length(varargin)>3
    fprintf('Too many input parameters! There are %d parameters but there should only be a maximum of 9', 6+length(varargin));
end


%% SelectedSets validation
%  First make sure there are valid sets for sorting

if iscell(EEGIn{1,1})
    HasMS = arrayfun(@(x) hasMicrostates(EEGIn{1,1}(x)), 1:numel(EEGIn)); % Call of subfunction to check that input datasets have microstate maps.
    HasDyn = arrayfun(@(x) isDynamicsSet(EEGIn{1,1}(x)), 1:numel(EEGIn)); % Call of subfunction to check if input datasets are dynamic sets.
    usingPublished = isPublishedSet(EEGIn{1,2}, EEGIn{1,2}.setname);
    AvailableSets = find(HasMS & ~HasDyn);
    HasChildren   = arrayfun(@(x) DoesItHaveChildren(EEGIn{1,1}(x)), AvailableSets);
    ChosenTemplate = EEGIn{1,2};
elseif isstruct(EEGIn{1,1})
    
    % Test the individual-level datasets to be sorted
    HasMS1 = cell2mat(arrayfun(@(x) clv_hasMicrostates(EEGIn{1,1}(x)), 1:numel(EEGIn{1,1}), 'UniformOutput',false)); % Call of subfunction to check that input datasets have microstate maps.
    HasDyn1 = cell2mat(arrayfun(@(x) clv_isDynamicsSet(EEGIn{1,1}(x)), 1:numel(EEGIn{1,1}), 'UniformOutput',false)); % Call of subfunction to check if input datasets are dynamic sets.
    usingPublished1 = cell2mat(arrayfun(@(x) clv_isPublishedSet(EEGIn{1,1}(x), {TemplateSet}), 1:numel(EEGIn{1,1}), 'UniformOutput',false));

    % Test the grand-mean template dataset.
    HasMS2 = cell2mat(arrayfun(@(x) clv_hasMicrostates(EEGIn{1,2}(x)), 1:numel(EEGIn{1,2}), 'UniformOutput',false));
    HasDyn2 = cell2mat(arrayfun(@(x) clv_isDynamicsSet(EEGIn{1,2}(x)), 1:numel(EEGIn{1,2}), 'UniformOutput',false));
    usingPublished2 = cell2mat(arrayfun(@(x) clv_isPublishedSet(EEGIn{1,2}(x), {TemplateSet}), 1:numel(EEGIn{1,2}), 'UniformOutput',false));

    AvailableSets_Ind = find(HasMS1 & ~HasDyn1 & ~usingPublished1);
    AvailableSets_Temp = find(HasMS2 & ~HasDyn2 & ~usingPublished2);

    HasChildren1 = cell2mat(arrayfun(@(x) clv_DoesItHaveChildren(EEGIn{1,1}(x)), AvailableSets_Ind, 'UniformOutput',false));
    HasChildren2 = cell2mat(arrayfun(@(x) clv_DoesItHaveChildren(EEGIn{1,2}(x)), AvailableSets_Temp, 'UniformOutput',false));
    ChosenTemplate = EEGIn{1,2}; 
    usingPublished = [sum(usingPublished1) sum(usingPublished2)];
end

AvailableIndSets = AvailableSets_Ind(~HasChildren1);  % Find the individual datasets - those that do not have children
AvailableMeanSets = AvailableSets_Temp(HasChildren2);  % Find the grandmean datasets - those with children.A

%% Verifying the input TemplateSet.

pluginPath = '/Users/bolger/Documents/Documents - MACPRO-BOLGER-01/MATLAB/eeglab2024.0/plugins/MICROSTATELAB2.1/Templates/';


if ischar(TemplateSet) && ~strcmp(TemplateSet, 'own')   % A published micro-state template set. Check that it exists.
    d = dir(pluginPath);
    if sum(contains({d.name}, TemplateSet))==1
        fprintf('The selected template set, %s, does exist in %s.\n', TemplateSet, pluginPath)
    elseif sum(contains({d.name}, TemplateSet))==0
        fprintf('The defined template set, %s, does not exist or is not available. \n', TemplateSet);
    end
end

% Need to add code to verify if a dataset is selected as a template set.

%% If datasets are being sorted based on a grand-mean dataset. Need to ensure that the datasets being sorted are children of the grand-mean dataset.

if strcmp(TemplateSet, 'GMdataset')
    gmchildren = {ChosenTemplate.msinfo.children}; % Find the children of the grand-mean dataset.
    
    if sum(cell2mat(arrayfun(@(x) ~isfield(EEGIn{1,1}(x).msinfo, 'children'), 1:numel(EEGIn{1,1}), 'UniformOutput',false))) == numel(EEGIn{1,1})

        arrayfun(@(x) fprintf('The current dataset to be sorted : %s is an individual-level dataset.\n', EEGIn{1,1}(x).setname), 1:numel(EEGIn{1,1}), 'UniformOutput', false)
        fprintf('Checking if is one of the children of the grand-mean template set...\n')
        ischildren = cell2mat(arrayfun(@(x) sum(ismember(gmchildren{1,1}, EEGIn{1,1}(x).setname)), 1:numel(EEGIn{1,1}), 'UniformOutput', false));  % Comparing the children of template dataset with datasets to sorted.
        
        if sum(ischildren)== numel(EEGIn{1,1})
            fprintf('The individual set, %s, is a child of the grand-mean template dataset.\n',EEGIn{1,1}.setname)
        end

    end

end

%% Verify the number of classes defined.
%  Check that the number of classes defined falls within the
%  msinfo.ClustPar.MinClasses and msinfo.ClustPar.MaxClasses
%  TemplateClasses cannot be set to "all" if TemplateSet is set "own".

if matches('own', TemplateSet) && matches('all', TemplateClassNum)
    temperror_Message = sprintf('The TemplateClassNum argument cannot be set to own if the TemplateSet argurment is all.\n '); 
    error(temperror_Message); 
end

%% Verifications prior to sorting and sorting of microstate maps.
%  If TemplateSet is defined as "own", verify that the template dataset
%  chosen is already sorted itself;

SortMode = ChosenTemplate.msinfo.MSMaps(TemplateClassNum).SortMode;
TemplateMinClasses = ChosenTemplate.msinfo.ClustPar.MinClasses;
TemplateMaxClasses = ChosenTemplate.msinfo.ClustPar.MaxClasses;

if matches('none', SortMode)
    errorMessage = sprintf('Selected template solution %i of template set %s is unsorted. Please sort before using as a template solution.', TemplateClasses, TemplateName);
    error(errorMessage);
end

if Stepwise

    for i=1:length(selectSets)
        fprintf('Sorting dataset %i of %i\n', i, numel(selectSets));
        sIndex = SelectSets(i);
        if iscell(EEGIn)
            MSMaps = StepwiseSort(EEGIn{1,1}(sIndex).msinfo.MSMaps, EEGIn{1,1}(sIndex).msinfo.ClustPar, EEGIn{1,1}(sIndex).setname, TemplateClasses, IgnorePolarity);
            if isempty(MSMaps); continue;   end
            EEGIn(sIndex).msinfo.MSMaps = MSMaps;
        elseif isstruct(EEGIn)
            MSMaps = StepwiseSort(EEGIn(sIndex).msinfo.MSMaps, EEGIn(sIndex).msinfo.ClustPar, EEGIn(sIndex).setname, TemplateClasses, IgnorePolarity);
            if isempty(MSMaps); continue;   end
            EEGIn(sIndex).msinfo.MSMaps = MSMaps;
        end
    end

else  % if not Stepwise

    for i = 1:length(selectSets)
        
        fprintf('Sorting dataset %i of %i\n', i, numel(selectSets));
        sIndex = selectSets(i);
        
        
        if strcmpi(ClassNum, 'all')

            if iscell(EEGIn)       
                sortClasses = EEGIn{1,1}(sIndex).msinfo.ClustPar.MinClasses:EEGIn{1,1}(sIndex).msinfo.ClustPar.MaxClasses;
            elseif isstruct(EEGIn)  
                 sortClasses = EEGIn(sIndex).msinfo.ClustPar.MinClasses:EEGIn(sIndex).msinfo.ClustPar.MaxClasses;
            end

        else
            sortClasses = ClassNum;
        end

        for nCnt = 1: numel(sortClasses)

            n = sortClasses(nCnt);

            % skip class number if the current set does not contain the
            % current cluster solution
            if (n > EEGIn{1,1}(sIndex).msinfo.ClustPar.MaxClasses) || (n < EEGIn{1,1}(sIndex).msinfo.ClustPar.MinClasses)
                fprintf('The current class number, %d, is not a member of the current cluster solution. \n', n);
            else
                fprintf('The current class number, %d, is a member of the current cluster solution. \n', n);
            end

            % find the number of template classes to use
            if ~strcmpi(TemplateClassNum, 'all')
                TemplateClassesToUse = TemplateClassNum;
            else
                if n < TemplateMinClasses
                    TemplateClassesToUse = TemplateMinClasses;
                elseif n > TemplateMaxClasses
                    TemplateClassesToUse = TemplateMaxClasses;
                else
                    TemplateClassesToUse = n;
                end
            end

            % skip class number if own maps are being sorted by the
            % same solution
            if strcmpi(TemplateSet, 'own') && (n == TemplateClassesToUse)
                fprintf('Skipping class number as own maps are being sorted by the same solution. \n')
            end

            % or if a mean set is being sorted by itself
            if ~strcmpi(TemplateSet, 'own') && sum(usingPublished)==0
                if strcmp(ChosenTemplate.setname, EEGIn{1,1}(sIndex).setname) && (n == TemplateClassesToUse)
                    fprintf('Skipping class number as the mean set is being sorted by itself. \n')
                end
            end

            if max(n, TemplateClassesToUse) >= 10 && (~license('test','optimization_toolbox') || isempty(which('intlinprog')))
                warning(['Sorting using 10 or more classes requires the Optimization toolbox. ' ...
                    'Please install the toolbox using the Add-On Explorer. Skipping large cluster solutions...']);
                
            end

            if ~strcmpi(TemplateSet, 'own')
                % compare number of channels in selected set and template set -
                % convert whichever set has more channels to the channel
                % locations of the other

                MapsToSort = zeros(1, n, min(EEGIn{1,1}(sIndex).nbchan, ChosenTemplate.nbchan));
                [LocalToGlobal, GlobalToLocal] = MakeResampleMatrices(EEGIn{1,1}(sIndex).chanlocs,ChosenTemplate.chanlocs);
                if EEGIn{1,1}(sIndex).nbchan > ChosenTemplate.nbchan
                    MapsToSort(1,:,:) = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps * LocalToGlobal';
                    TemplateMaps      = ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).Maps;
                else
                    MapsToSort(1,:,:) = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps;
                    TemplateMaps      = ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).Maps * GlobalToLocal';
                end
            else
                MapsToSort = zeros(1, n, EEGIn{1,1}(sIndex).nbchan);
                MapsToSort(1,:,:) = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps;
                ChosenTemplate    = EEGIn{1,1}(sIndex);
                TemplateMaps      = EEGIn{1,1}(sIndex).msinfo.MSMaps(TemplateClassesToUse).Maps;
            end

            % If the template set has unassigned maps, remove them (only
            % base sorting on assigned maps)
            nAssignedLabels = sum(~arrayfun(@(x) all(ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).ColorMap(x,:) == [.75 .75 .75]), 1:TemplateClassesToUse));
            if nAssignedLabels < TemplateClassesToUse && nAssignedLabels > 0
                TemplateMaps(nAssignedLabels+1:end,:) = [];
            end

            % Sort
            [~,SortOrder, SpatialCorrelation, polarity] = ArrangeMapsBasedOnMean(MapsToSort,TemplateMaps);
            EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps(SortOrder(SortOrder <= n),:);
            EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Maps .* repmat(polarity',1,numel(EEGIn{1,1}(sIndex).chanlocs));

            % Update map labels and colors
            [Labels,Colors] = UpdateMicrostateLabels(EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Labels,ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).Labels,SortOrder,ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).ColorMap);
            EEGIn{1,1}(sIndex).msinfo.MSMaps(n).Labels = Labels;
            EEGIn{1,1}(sIndex).msinfo.MSMaps(n).ColorMap = Colors;

            % Update individual explained variance order
            if numel(EEGIn{1,1}(sIndex).msinfo.MSMaps(n).ExpVar) > 1
                EEGIn{1,1}(sIndex).msinfo.MSMaps(n).ExpVar = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).ExpVar(SortOrder(SortOrder <= n));
            end

            % Update shared variance order
            if isfield(EEGIn{1,1}(sIndex).msinfo.MSMaps(n), 'SharedVar')
                EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SharedVar = EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SharedVar(SortOrder(SortOrder <= n));
            end

            if strcmpi(TemplateSet, 'own')
                EEGIn(sIndex).msinfo.MSMaps(n).SortMode = 'own template maps';
                EEGIn(sIndex).msinfo.MSMaps(n).SortedBy = sprintf('%s->%s (%i classes)', ...
                    EEGIn(sIndex).msinfo.MSMaps(TemplateClassesToUse).SortedBy, EEGIn(sIndex).setname, TemplateClassesToUse);
            else
                if usingPublished
                    EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SortMode = 'published template maps';
                    EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SortedBy = sprintf('%s (%i classes)', ChosenTemplate.setname, TemplateClassesToUse);
                else
                    if strcmp(TemplateSet, EEGIn{1,1}(sIndex).setname)
                        EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SortMode = 'own template maps';
                    else
                        EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SortMode = 'mean template maps';
                    end
                    EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SortedBy = sprintf('%s->%s (%i classes)', ...
                        ChosenTemplate.msinfo.MSMaps(TemplateClassesToUse).SortedBy, ChosenTemplate.setname, TemplateClassesToUse);
                end
            end

            EEGIn{1,1}(sIndex).msinfo.MSMaps(n).SpatialCorrelation = SpatialCorrelation;
        end

        EEGIn{1,1}(sIndex).saved = 'no';

        savepath = EEGIn{1,1}(sIndex).filepath;
        EEG = pop_saveset(EEGIn{1,1}(sIndex), 'filename', [EEGIn{1,1}(sIndex).setname,'_sortedPublished'], 'filepath', savepath);
      
    end
end

EEGout = EEGIn{1,1}(selectSets);
CurrentSet = selectSets;

% Save the datasets with the sorted MS maps. 


end % end of function

%% ********************** Sub-Functions **********************************

function hasDyn = isDynamicsSet(in)
hasDyn = false;
% check if set includes msinfo
if ~isfield(in,'msinfo')
    return;
end
% check if set has FitPar
if ~isfield(in.msinfo, 'FitPar')
    return;
end
% check if FitPar contains Rectify/Normalize parameters
if ~isfield(in.msinfo.FitPar, 'Rectify')
    return;
else
    hasDyn = true;
end
end

%% ***********************************************************

function hasMS = hasMicrostates(in)
hasMS = false;

% check if set includes msinfo
if ~isfield(in,'msinfo')
    return;
end

% check if msinfo is empty
if isempty(in.msinfo)
    return;
else
    hasMS = true;
end
end

%% ***************************************************************
function Answer = DoesItHaveChildren(in)
Answer = false;
if ~isfield(in,'msinfo')
    return;
end

if ~isfield(in.msinfo,'children')
    return
else
    Answer = true;
end
end

%% *****************************************************************

function isPublished = isPublishedSet(in, templateNames)
isPublished = false;
if isempty(in.setname)
    return;
end

if matches(in.setname, templateNames)
    isPublished = true;
end
end

%% *****************************************************************

function containsChild = checkSetForChild(AllEEGIn, SetsToSearch, childSetName)
containsChild = false;
if isempty(SetsToSearch)
    return;
end

% find which sets have children
HasChildren = arrayfun(@(x) isfield(AllEEGIn(x).msinfo, 'children'), SetsToSearch);
% if none of the sets to search have children, the child set could not
% be found
if ~any(HasChildren)
    return;
end

% find non-empty children fields
nonempty = arrayfun(@(x) ~isempty(AllEEGIn(x).msinfo.children), SetsToSearch(HasChildren));
if ~any(nonempty)
    return;
end
HasChildren = HasChildren(nonempty);

% search the children of all the mean sets for the child set name
containsChild = any(arrayfun(@(x) matches(childSetName, AllEEG(x).msinfo.children), SetsToSearch(HasChildren)));

% if the child cannot be found, search the children of the children
if ~containsChild
    setnames = {AllEEGIn.setname};
    isEmpty = cellfun(@isempty,setnames);
    if any(isEmpty)
        setnames(isEmpty) = {''};
    end
    childSetIndices = unique(cell2mat(arrayfun(@(x) find(matches(setnames, AllEEGIn(x).msinfo.children)), SetsToSearch(HasChildren), 'UniformOutput', false)));
    containsChild = checkSetForChild(AllEEG, childSetIndices, childSetName);
end

end











