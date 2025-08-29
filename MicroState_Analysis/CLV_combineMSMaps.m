
function EEGout = CLV_combineMSMaps(EEGIn, selectSets, MeanName, IgnorePolarity, ShowMaps, msfigsave)
%% Combine micro-state maps across participant resting-state blocks, groups and conditions
%  Based on the MICROSTATELAB function pop_CombMSMaps().
%
%  This function carries out the identification of the mean and grand mean
%  micro-state maps. The individual-level micro-state maps are re-ordered
%  such that the shared variance across participants is maximized.
%
%   Input parameters: 
%   EEGIn : the datasets with individual-level microstate information
%   (structure)
%   selectSets : the indices of those datasets to analyse (vector)
%   MeanName : the name to give to the output dataset with the averaged and
%   sorted microstate map information. (character array)
%   IgnorePolarity : true or false if polarity was ignored when calculating
%   the microstate maps. (0 or 1, default is 1)
%   ShowMaps : true or false if to create figure of the mean and sorted
%   microstate maps. (0 or 1, default is 1)
%   msfigsave : the path in which to save the mean microstate figures.
%   (character array)
%
%   Output parameters :
%   EEGout : dataset with mean and sorted microstate maps for the current
%   condition or group. (structure)
%
% *****************************************************************************************

%% Get the microstate parameters from the msinfo field of the selected datasets.
%  Verify the microstate parameters and get channel labels.

MinClasses     = EEGIn(selectSets(1)).msinfo.ClustPar.MinClasses;       fprintf("The minimum number of MS classes is %d.\n", MinClasses)
MaxClasses     = EEGIn(selectSets(1)).msinfo.ClustPar.MaxClasses;       fprintf("The maximum number of classes is %d.\n", MaxClasses)
tempIPolarity  = EEGIn(selectSets(1)).msinfo.ClustPar.IgnorePolarity;   fprintf("The choice for ignore polarity is %d.\n", tempIPolarity)
GFPPeaks       = EEGIn(selectSets(1)).msinfo.ClustPar.GFPPeaks;         fprintf("The choice for applying GFP peaks is %d.\n", GFPPeaks)


if ~isfield(EEGIn(selectSets(1)).msinfo.ClustPar,'UseEMD')
    UseEMD = false;
else
    UseEMD = EEGIn(selectSets(1)).msinfo.ClustPar.UseEMD;
end

allchans  = { };
children  = cell(length(selectSets),1);
keepindex = 0;

for index = 1:length(selectSets)
    if  MinClasses     ~= EEGIn(selectSets(index)).msinfo.ClustPar.MinClasses || ...
            MaxClasses     ~= EEGIn(selectSets(index)).msinfo.ClustPar.MaxClasses || ...
            tempIPolarity  ~= EEGIn(selectSets(index)).msinfo.ClustPar.IgnorePolarity || ...
            GFPPeaks       ~= EEGIn(selectSets(index)).msinfo.ClustPar.GFPPeaks
        errordlg2('Microstate parameters differ between datasets','Identify mean microstate maps error');
        return;
    end

    children(index) = {EEGIn(selectSets(index)).setname};
    tmpchanlocs = EEGIn(selectSets(index)).chanlocs;
    tmpchans = { tmpchanlocs.labels };
    allchans = unique_bc([ allchans {tmpchanlocs.labels}]);

    if length(allchans) == length(tmpchans)
        keepindex = index;
    end
end

if keepindex
    tmpchanlocs = EEGIn(selectSets(keepindex)).chanlocs;
end

msinfo.children = children;
msinfo.ClustPar = EEGIn(selectSets(1)).msinfo.ClustPar;

%% Create the combined maps.

for n = MinClasses:MaxClasses

    MapsToSort = nan(numel(selectSets),n,numel(tmpchanlocs));
    % Here we go to the common set of channels
    for index = 1:length(selectSets)
        LocalToGlobal = MakeResampleMatrices(EEGIn(selectSets(index)).chanlocs,tmpchanlocs);
        MapsToSort(index,:,:) = L2NormDim(EEGIn(selectSets(index)).msinfo.MSMaps(n).Maps * LocalToGlobal',2);
    end
    % We sort out the stuff
    [BestMeanMap,~,ExpVar,SharedVar] = PermutedMeanMaps(MapsToSort,~IgnorePolarity,tmpchanlocs,[],UseEMD); % Call of function PermutedMeanMaps

    msinfo.MSMaps(n).Maps = BestMeanMap;
    msinfo.MSMaps(n).ExpVar = ExpVar;
    msinfo.MSMaps(n).SharedVar = SharedVar;
    msinfo.MSMaps(n).ColorMap = repmat([.75 .75 .75], n, 1);
    for j = 1:n
        msinfo.MSMaps(n).Labels{j} = sprintf('MS_%i.%i', n,j);
    end
    msinfo.MSMaps(n).SortMode = 'none';
    msinfo.MSMaps(n).SortedBy = 'none';
    msinfo.MSMaps(n).SpatialCorrelation = [];

end

EEGout = eeg_emptyset();
EEGout.chanlocs = tmpchanlocs;
EEGout.data = zeros(numel(EEGout.chanlocs),MaxClasses,MaxClasses);
EEGout.msinfo = msinfo;

for n = MinClasses:MaxClasses
    EEGout.data(:,1:n,n) = msinfo.MSMaps(n).Maps';
end

EEGout.setname     = MeanName;
EEGout.nbchan      = size(EEGout.data,1);
EEGout.trials      = size(EEGout.data,3);
EEGout.pnts        = size(EEGout.data,2);
EEGout.srate       = 1;
EEGout.xmin        = 1;
EEGout.times       = 1:EEGout.pnts;
EEGout.xmax        = EEGout.times(end);

%% 

if ShowMaps
    fig_h = arrayfun(@(a) pop_ShowIndMSMaps(a, 1:numel(a), 'Classes', MinClasses:MaxClasses, 'Visible', 0), EEGout, 'UniformOutput', false);
    
    currfigName = [MeanName,'_microstate_maps.png'];
    
    fig_fullPath = fullfile(msfigsave, currfigName);
    saveas(fig_h{1,1}, fig_fullPath);
end

end