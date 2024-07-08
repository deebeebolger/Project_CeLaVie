%% Programmer: D. Bolger                 Date: 22-05-2024
%  Function to run the AMICA algorithm. 
%  AMICA version: 1.7
%  Inputs:
%  - EEG: Input EEG structure 
%  - RELAX_cfg: Configuration structure for RELAX pipeline.
%  Outputs:
%  - Outeeg: Output EEG structure 
%  - wICA:   wavelet thresholded ICs
%  - A:      demixing matrix.
%  - IC:     non-wavelet ICs.
%**************************************************************************

function [Outeeg, wICA, A, W, IC] = CLV_runamica(EEG, RELAX_cfg)

% Maybe run this for the electrodes of type EEG.
iseeg = ismember({EEG.chanlocs.type}, 'EEG');


% Define amica output directory.
amica_outdir = fullfile(RELAX_cfg.myPath,'amicaouttmp', filesep);
mkdir(amica_outdir);

% Get the rank of the data. 
tmprank = getrank(EEG.data(:,1:min(3000, size(EEG.data,2))));

% Define the parameters.
numprocs    = 1;      % The number of nodes
max_threads = 24;
num_models  = 1;      % The number of models of mixture ICA.
max_iter    = 1000;   % The max number of learning steps.

% Run amica
[weights, sphere, mods] = runamica15(EEG.data, 'num_models', num_models, 'outdir', amica_outdir, ...
    'numprocs', numprocs, 'max_threads', max_threads, 'max_iter', max_iter);

end 