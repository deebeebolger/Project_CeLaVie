function EEGOut = CLV_findMS_maps(DataIn, varargin)
%%******************************************************************
% Date : June 2025        Programmed by : D. Bolger
% Function to calculate the microstate maps at the individual,
% subject-level. This script is based on the findMS_maps() scripts from the
% Microstatelab toolbox for EEGLAB (Koenig & Aryan, 2023). 
% Reference : "MICROSTATELAB: The EEGLAB toolbox for resting-state microstate 
% analysis by Thomas Koenig and Delara Aryan"
% It applies the modified kmeans algorithm to carry out Microstate
% clustering. This modified kmeans ensures that maps of opposite polarity
% are considered equivalent. 
% 


%% Parse inputs and perform initial validation
SelectedSets = varargin{1,1};
ClustPar = varargin{1,2};
TTFrD = varargin{1,3};

ShowMaps = true;
FailedSets = [];

%% Start MS segmentation on datasets

for counter = 1 : numel(SelectedSets)

    currSetIndx = SelectedSets(counter);
    fprintf('Carrying out microstate clustering on data %d : %s. \n', SelectedSets(counter), DataIn(currSetIndx).setname);

    % Distribute random sampling across temporal segments of data.
    nSegs = DataIn(currSetIndx).trials;   % In spontaneous data this is 1.
    if ~isinf(ClustPar.MaxMaps)
        MapsPerSegment = hist(ceil(double(nSegs) * rand(ClustPar.MaxMaps,1)),nSegs);
    else
        MapsPerSegment = inf(nSegs,1);  % An infinite number of maps per segment.
    end

    %% If DoTTFrD is false. Compute the maps to use for MS segmentation;

    MapsToUse = [];
    for s = 1:nSegs
        if ClustPar.GFPPeaks == 1
            gfp = std(DataIn(currSetIndx).data(:,:,s),1,1);                                                   % Calculate GFP for current dataset
            IsGFPPeak = find([false (gfp(1,1:end-2) < gfp(1,2:end-1) & gfp(1,2:end-1) > gfp(1,3:end)) false]);  % Find GFP peaks
            if numel(IsGFPPeak) > MapsPerSegment(s) && MapsPerSegment(s) > 0
                idx = randperm(numel(IsGFPPeak));
                IsGFPPeak = IsGFPPeak(idx(1:MapsPerSegment(s)));
            end
            MapsToUse = [MapsToUse DataIn(currSetIndx).data(:,IsGFPPeak,s)];
        else
            if (size(DataIn(currSetIndx).data,2) > ClustPar.MaxMaps) && MapsPerSegment(s) > 0
                idx = randperm(size(DataIn(currSetIndx).data,2));
                MapsToUse = [MapsToUse DataIn{1,currSetIndx}.data(:,idx(1:MapsPerSegment(s)),s)];
            else
                MapsToUse = [MapsToUse DataIn{1,currSetIndx}.data(:,:,s)];
            end
        end
    end

    if size(MapsToUse,2) < ClustPar.MaxClasses
        warning('Not enough data to cluster in set %s',DataIn{1,currSetIndx}.setname);
        FailedSets = [FailedSets,currSetIndx];
        continue;
    end

    %% Define flags and carry out MS clustering.

    flags = '';
    if ClustPar.IgnorePolarity == false
        flags = [flags 'p'];
    end
    if ClustPar.Normalize == true
        flags = [flags 'n'];
    end

    if isfield(ClustPar, 'UseEMD')
        if ClustPar.UseEMD == true
            flags = [flags 'e'];
        end
    end

    if ClustPar.UseAAHC == false
        for nClusters = ClustPar.MinClasses:ClustPar.MaxClasses
            [b_model,~,~,exp_var] = eeg_kMeans(MapsToUse',nClusters,ClustPar.Restarts,[],flags,DataIn(currSetIndx).chanlocs);

            msinfo.MSMaps(nClusters).Maps = double(b_model);
            msinfo.MSMaps(nClusters).ExpVar = double(exp_var);
            msinfo.MSMaps(nClusters).ColorMap = repmat([.75 .75 .75], nClusters, 1);
            for j = 1:nClusters
                msinfo.MSMaps(nClusters).Labels{j} = sprintf('MS_%i.%i',nClusters,j);
            end
            msinfo.MSMaps(nClusters).SortMode = 'none';
            msinfo.MSMaps(nClusters).SortedBy = '';
            msinfo.MSMaps(nClusters).SpatialCorrelation= [];
        end
    else
        [b_model,exp_var] = eeg_computeAAHC(double(MapsToUse'),ClustPar.MinClasses:ClustPar.MaxClasses,false, ClustPar.IgnorePolarity,ClustPar.Normalize);

        for nClusters = ClustPar.MinClasses:ClustPar.MaxClasses
            msinfo.MSMaps(nClusters).Maps = double(b_model{nClusters-ClustPar.MinClasses+1});
            msinfo.MSMaps(nClusters).ExpVar = double(exp_var{nClusters-ClustPar.MinClasses+1});
            msinfo.MSMaps(nClusters).ColorMap = repmat([.75 .75 .75], nClusters, 1);
            for j = 1:nClusters
                msinfo.MSMaps(nClusters).Labels{j} = sprintf('MS_%i.%i',nClusters,j);
            end
            msinfo.MSMaps(nClusters).SortMode = 'none';
            msinfo.MSMaps(nClusters).SortedBy = '';
            msinfo.MSMaps(nClusters).SpatialCorrelation= [];
        end
    end

    msinfo.ClustPar = ClustPar;
    DataIn(currSetIndx).msinfo = msinfo;
    DataIn(currSetIndx).saved = 'no';



end

% Remove sets that were not clustered
SelectedSets(FailedSets) = [];

%% Show maps
EEGOut = DataIn(SelectedSets);
nameSplit = cellfun(@(b) split(b,"_"),{EEGOut.setname}, 'UniformOutput',false);
savepath = fullfile(cd, 'Data', 'Microstate', 'Figures');

if ShowMaps
   
    fig_h = arrayfun(@(a) pop_ShowIndMSMaps(a, 1:numel(a), 'Classes', ClustPar.MinClasses:ClustPar.MaxClasses, 'Visible', 0), EEGOut, 'UniformOutput', false);
    
    for i = 1:numel(nameSplit)
        rsname = nameSplit{1,i}(find(contains(nameSplit{1,i},'restingstate')));
        if numel(rsname)>1
            rsname = rsname{2,1};
        end
        currfigName = [char(nameSplit{1,i}(1)),'-',rsname,'-microstate_maps.png'];
        fig_fullPath = fullfile(savepath, currfigName);
        saveas(fig_h{1,i},fig_fullPath);
    end
end



end
