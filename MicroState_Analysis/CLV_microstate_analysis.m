function CLV_microstate_analysis()
%% ***************************************************************************************
% Date : Juin 2025    Programmed by : D. BOLGER
% Script to carry out the microstate analysis of the CeLaVie data using the
% Microstatelab toolbox, a EEGLAB plugin (Kalburgi et al, 2024). 
% The script carries out the following steps in the analysis of microstates
% of spontaneous data. 
%
%
%% Load in the datasets

fnameAll = {};
fpathAll = {};

[filename, filepath] = uigetfile({'*.bdf; *.set' 'BDF or eeglab set file';'*.bdf' 'BDF'; '*.set' 'eeglab set file'}, 'Choose a BDF or .set data file ', 'multiselect', 'on');
[parentdir, ~,~] = fileparts(filepath);
[~,~,ext] = fileparts(filename);
main_dir = fileparts(parentdir);
fnameAll = [fnameAll; {filename}];
fpathAll = [fpathAll; {filepath}];
fpathAll = repmat(fpathAll,1,size(filename,2));   % Add repeated copies of filepath to have same dimensions as filenames.


dataIn = cellfun(@(x,y) pop_loadset(x,y), fnameAll{1,1}, fpathAll, 'UniformOutput',false); % Load in the datasets. 



%% Carry out the microstate segmentation of the individual datasets.

Clustpar = [];                  % Cluster paramters to define and call in the pop_FindMSMaps () function.
Clustpar.UseAAHC = false;       % Carries out kmeans, if true it carries out AAHC
Clustpar.MinClasses = 4;        % Minimum number of clusters to identify
ClustPar.MaxClasses = 7;        % Maximum number of clusters to identify
ClustPar.Restarts 20;           % Number of times kmeans algorithm is restarted (this is ignored if using AAHC)
ClustPar.MaxMaps = inf;         % Maximum number of data samples to use to identify clusters.
ClustPar.GFPPeaks = true;       % Whether clustering should be limited to GFP peaks.
ClustPar.IgnorePolarity = true; % whether maps with inverted polarity should be considered part of the same cluster.
ClustPar.Normalize = true;      % set to false if using AAHC.

pop_FindMSMaps(dataIn, 1, 'ClustPar', ClustPar)





end