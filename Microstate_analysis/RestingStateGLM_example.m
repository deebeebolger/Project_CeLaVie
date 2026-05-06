
%% Example: Complete Resting-State EEG Analysis

clear; clc;

fprintf('========================================\n');
fprintf('  Resting-State EEG Spatial GLM\n');
fprintf('========================================\n\n');

%% Step 1: Load or generate resting-state EEG data

% Option A: Load real data
% EEG = pop_loadset('resting_state.set');

% Option B: Generate synthetic resting-state data
EEG = generate_synthetic_resting_eeg();

%% Step 2: Extract microstate maps

% Option 1: Use predefined canonical microstates
n_microstates = 4;
microstate_maps = create_canonical_microstates(EEG.chanlocs, n_microstates);

% Option 2: Extract from the data using modified K-means
% microstate_maps = extract_microstates_kmeans(EEG, n_microstates);

%% Step 3: Initialize Spatial GLM

glm = RestingStateGLM(EEG, microstate_maps, EEG.chanlocs, EEG.srate);

%% Step 4: Fit GLM using different methods

fprintf('\n--- Fitting GLM ---\n');

% Method 1: Standard OLS (fast but noisy)
% glm.fit_continuous_glm('ols');

% Method 2: Tapered sliding window (smooth, recommended)
params.window_size = round(0.05 * EEG.srate); % 50 ms window
glm.fit_continuous_glm('tapered_sliding', params);

% Method 3: AAHC (discrete microstate assignment)
% params.smoothing = 5;
% glm.fit_continuous_glm('aahc', params);

%% Step 5: Extract microstate parameters

glm.extract_microstate_parameters(10); % Minimum 10 ms duration

%% Step 6: Spectral analysis

freq_bands.delta = [1 4];
freq_bands.theta = [4 8];
freq_bands.alpha = [8 13];
freq_bands.beta = [13 30];
glm.analyze_spectral_properties(freq_bands);

%% Step 7: Spatial autocorrelation analysis

glm.create_spatial_weights('knn', 4);
glm.test_spatial_autocorrelation();

%% Step 8: Stationarity analysis

glm.analyze_stationarity(10); % 10-second segments

%% Step 9: Visualizations

% Coefficient time courses
glm.plot_coefficient_timecourses([0, 10]);

% Microstate topographies
glm.plot_microstate_topographies();

% Transition network
glm.plot_tran