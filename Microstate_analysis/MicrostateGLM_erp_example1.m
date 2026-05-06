%% Example: Microstate GLM for Event-Related Potentials
%  Complete example of microstateGLM for ERPs

clear; clc;

% Check for EEGLAB
if ~exist('eeglab', 'file')
    warning('EEGLAB not found. Using synthetic data.');
    USE_REAL_DATA = false;
else
    USE_REAL_DATA = true;
end

%% Step 1: Load or create EEG data
if USE_REAL_DATA
    % Load real EEG data (replace with your file)
    % EEG = pop_loadset('your_erp_data.set');
    
    % For demonstration, create synthetic ERP data
    EEG = create_synthetic_erp_data();
else
    EEG = create_synthetic_erp_data();
end

%% Step 2: Define microstate maps
% Option A: Use predefined canonical microstates
n_microstates = 4;
microstate_maps = create_canonical_microstates(EEG.chanlocs, n_microstates);

% Option B: Extract microstates from resting-state data
% microstate_maps = extract_microstates_from_data(resting_EEG, 4);

%% Step 3: Initialize and fit GLM
fprintf('\n========================================\n');
fprintf('   Microstate GLM Analysis\n');
fprintf('========================================\n\n');
fprintf('Data dimensions:\n');
fprintf('  Channels: %d\n', EEG.nbchan);
fprintf('  Time points: %d\n', EEG.pnts);
fprintf('  Trials: %d\n', EEG.trials);
fprintf('  Microstates: %d\n\n', n_microstates);

% Create GLM object
glm = MicrostateGLM(EEG, microstate_maps);

% Fit using OLS
glm.fit_glm('ols');

% Optional: Apply temporal smoothing
glm.smooth_coefficients(10, 'gaussian');

%% Step 4: Compute statistics
times = EEG.times; % Time vector in ms
baseline_window = [min(times), 0]; % Pre-stimulus baseline
stats = glm.compute_statistics(baseline_window, times);

%% Step 5: Visualize results

% Plot coefficient time courses
glm.plot_coefficient_timecourses(times, false);

% Plot topographies at key time points
key_timepoints = [0, 100, 200, 300, 400]; % ms
glm.plot_topographies(key_timepoints, times);

% Additional visualization: Heatmap of all coefficients
figure('Position', [100, 100, 1000, 400]);
subplot(1, 2, 1);
if glm.n_trials == 1
    beta_avg = glm.beta_timecourse;
else
    beta_avg = mean(glm.beta_timecourse, 3);
end
imagesc(times, 1:n_microstates, beta_avg);
xlabel('Time (ms)');
ylabel('Microstate');
title('Coefficient Time Courses (Heatmap)');
colorbar;
colormap(jet);

subplot(1, 2, 2);
imagesc(times, 1:n_microstates, stats.beta_zscore);
xlabel('Time (ms)');
ylabel('Microstate');
title('Z-scored Coefficients');
colorbar;
colormap(jet);
hold on;
plot([0 0], ylim, 'w--', 'LineWidth', 2); % Stimulus onset

%% Step 6: Statistical testing (optional)
% Perform cluster-based permutation test
fprintf('\nPerforming cluster-based permutation test...\n');
perm_results = cluster_permutation_test(glm, baseline_window, times);

%% Step 7: Export results
glm.export_results('microstate_glm_results.mat');

fprintf('\nAnalysis complete!\n');