%% GLM for Microstate Coefficient Time Courses
% Solves: EEG_data(t) = Σ β_i(t) * Microstate_i + ε(t)
% where β_i(t) are time-varying coefficients

%% Core GLM solver

classdef MicrostateGLM < handle
    properties
        EEG                 % EEGLAB structure or EEG data
        microstate_maps     % [n_channels x n_microstates] spatial templates
        n_channels          % Number of EEG channels
        n_microstates       % Number of microstate classes
        n_timepoints        % Number of time points
        n_trials            % Number of trials/epochs
        beta_timecourse     % [n_microstates x n_timepoints] or [n_microstates x n_timepoints x n_trials]
        fitted_data         % Reconstructed EEG data
        residuals           % Residual variance
        GEV                 % Global Explained Variance
        R2_timecourse       % R² at each time point
        chanlocs            % Channel locations
    end
    
    methods
        function obj = MicrostateGLM(EEG, microstate_maps)
            % Constructor
            % EEG: EEGLAB structure or [n_channels x n_timepoints x n_trials]
            % microstate_maps: [n_channels x n_microstates]
            
            if isstruct(EEG)
                % EEGLAB structure
                obj.EEG = EEG.data;
                obj.chanlocs = EEG.chanlocs;
                obj.n_channels = EEG.nbchan;
                obj.n_timepoints = EEG.pnts;
                obj.n_trials = EEG.trials;
            else
                % Raw data matrix
                obj.EEG = EEG;
                [obj.n_channels, obj.n_timepoints, obj.n_trials] = size(EEG);
                obj.chanlocs = [];
            end
            
            obj.microstate_maps = microstate_maps;
            obj.n_microstates = size(microstate_maps, 2);
            
            % Normalize microstate maps (L2 norm)
            for i = 1:obj.n_microstates
                obj.microstate_maps(:, i) = obj.microstate_maps(:, i) / ...
                    norm(obj.microstate_maps(:, i));
            end
        end
        
        function fit_glm(obj, method, regularization)
            % Fit GLM to obtain coefficient time courses
            % method: 'ols', 'ridge', 'lasso', 'spatial_ridge'
            % regularization: lambda parameter (for ridge/lasso)
            
            if nargin < 2, method = 'ols'; end
            if nargin < 3, regularization = 0; end
            
            fprintf('Fitting GLM using %s method...\n', upper(method));
            
            % Initialize coefficient matrix
            if obj.n_trials == 1
                obj.beta_timecourse = zeros(obj.n_microstates, obj.n_timepoints);
                obj.fitted_data = zeros(obj.n_channels, obj.n_timepoints);
                obj.residuals = zeros(obj.n_channels, obj.n_timepoints);
            else
                obj.beta_timecourse = zeros(obj.n_microstates, obj.n_timepoints, obj.n_trials);
                obj.fitted_data = zeros(obj.n_channels, obj.n_timepoints, obj.n_trials);
                obj.residuals = zeros(obj.n_channels, obj.n_timepoints, obj.n_trials);
            end
            
            obj.R2_timecourse = zeros(obj.n_timepoints, obj.n_trials);
            
            A = obj.microstate_maps; % Design matrix [n_channels x n_microstates]
            
            % Solve for each time point and trial
            for trial = 1:obj.n_trials
                for t = 1:obj.n_timepoints
                    % Extract data at time t
                    if obj.n_trials == 1
                        y = obj.EEG(:, t);
                    else
                        y = obj.EEG(:, t, trial);
                    end
                    
                    % Solve GLM based on method
                    switch lower(method)
                        case 'ols'
                            % Ordinary Least Squares: β = (A'A)^(-1) A'y
                            beta = (A' * A) \ (A' * y);
                            
                        case 'ridge'
                            % Ridge regression: β = (A'A + λI)^(-1) A'y
                            lambda = regularization;
                            beta = (A' * A + lambda * eye(obj.n_microstates)) \ (A' * y);
                            
                        case 'lasso'
                            % Lasso (L1 regularization)
                            lambda = regularization;
                            beta = lasso_solver(A, y, lambda);
                            
                        case 'nnls'
                            % Non-negative least squares
                            beta = lsqnonneg(A, y);
                            
                        case 'robust'
                            % Robust regression (iteratively reweighted LS)
                            beta = robustfit(A, y, 'huber', [], 'off');
                            
                        case 'spatial_ridge'
                            % Spatially regularized ridge
                            beta = obj.spatial_ridge_solver(A, y, regularization);
                    end
                    
                    % Store coefficients
                    if obj.n_trials == 1
                        obj.beta_timecourse(:, t) = beta;
                    else
                        obj.beta_timecourse(:, t, trial) = beta;
                    end
                    
                    % Compute fitted values and residuals
                    y_fitted = A * beta;
                    residual = y - y_fitted;
                    
                    if obj.n_trials == 1
                        obj.fitted_data(:, t) = y_fitted;
                        obj.residuals(:, t) = residual;
                    else
                        obj.fitted_data(:, t, trial) = y_fitted;
                        obj.residuals(:, t, trial) = residual;
                    end
                    
                    % Compute R²
                    SS_res = sum(residual.^2);
                    SS_tot = sum((y - mean(y)).^2);
                    obj.R2_timecourse(t, trial) = 1 - SS_res / SS_tot;
                end
                
                if mod(trial, 10) == 0
                    fprintf('  Processed trial %d/%d\n', trial, obj.n_trials);
                end
            end
            
            % Compute Global Explained Variance (GEV)
            obj.compute_GEV();
            
            fprintf('GLM fitting complete!\n');
            fprintf('Global Explained Variance (GEV): %.2f%%\n', obj.GEV * 100);
        end
        
        function beta = spatial_ridge_solver(obj, A, y, lambda)
            % Spatially-regularized ridge regression
            % Penalizes spatial non-smoothness of coefficients
            
            % Create spatial penalty matrix based on microstate similarity
            D = pdist2(A', A');  % Dissimilarity between microstates
            L = diag(sum(D, 2)) - D;  % Graph Laplacian
            
            % Solve: β = (A'A + λL)^(-1) A'y
            beta = (A' * A + lambda * L) \ (A' * y);
        end
        
        function compute_GEV(obj)
            % Compute Global Explained Variance
            
            if obj.n_trials == 1
                variance_orig = sum(var(obj.EEG, 0, 2));
                variance_res = sum(var(obj.residuals, 0, 2));
            else
                % Average across trials
                variance_orig = sum(var(mean(obj.EEG, 3), 0, 2));
                variance_res = sum(var(mean(obj.residuals, 3), 0, 2));
            end
            
            obj.GEV = (variance_orig - variance_res) / variance_orig;
        end
        
        function smooth_coefficients(obj, window_size, method)
            % Temporal smoothing of coefficient time courses
            % window_size: smoothing window in samples
            % method: 'moving_average', 'gaussian', 'savitzky_golay'
            
            if nargin < 3, method = 'moving_average'; end
            
            fprintf('Smoothing coefficient time courses (%s, window=%d)...\n', ...
                    method, window_size);
            
            for ms = 1:obj.n_microstates
                for trial = 1:obj.n_trials
                    if obj.n_trials == 1
                        signal = obj.beta_timecourse(ms, :);
                    else
                        signal = squeeze(obj.beta_timecourse(ms, :, trial));
                    end
                    
                    switch lower(method)
                        case 'moving_average'
                            smoothed = movmean(signal, window_size);
                            
                        case 'gaussian'
                            sigma = window_size / 4;
                            smoothed = imgaussfilt(signal, sigma);
                            
                        case 'savitzky_golay'
                            order = 3;
                            smoothed = sgolayfilt(signal, order, window_size);
                            
                        case 'lowpass'
                            % Low-pass filter
                            fc = 1 / window_size; % Cutoff frequency
                            [b, a] = butter(4, fc);
                            smoothed = filtfilt(b, a, signal);
                    end
                    
                    if obj.n_trials == 1
                        obj.beta_timecourse(ms, :) = smoothed;
                    else
                        obj.beta_timecourse(ms, :, trial) = smoothed;
                    end
                end
            end
        end
        
        function stats = compute_statistics(obj, baseline_window, times)
            % Compute statistics for coefficient time courses
            % baseline_window: [start, end] in ms for baseline correction
            % times: time vector in ms
            
            if nargin < 3
                times = 1:obj.n_timepoints;
            end
            
            % Find baseline indices
            if nargin > 1 && ~isempty(baseline_window)
                baseline_idx = times >= baseline_window(1) & times <= baseline_window(2);
            else
                baseline_idx = 1:round(obj.n_timepoints * 0.1); % First 10%
            end
            
            % Initialize statistics structure
            stats.mean_beta = zeros(obj.n_microstates, obj.n_timepoints);
            stats.sem_beta = zeros(obj.n_microstates, obj.n_timepoints);
            stats.baseline_mean = zeros(obj.n_microstates, 1);
            stats.baseline_std = zeros(obj.n_microstates, 1);
            stats.beta_zscore = zeros(obj.n_microstates, obj.n_timepoints);
            
            for ms = 1:obj.n_microstates
                if obj.n_trials == 1
                    beta_trace = obj.beta_timecourse(ms, :);
                    stats.mean_beta(ms, :) = beta_trace;
                    stats.sem_beta(ms, :) = zeros(1, obj.n_timepoints);
                else
                    beta_trace = squeeze(obj.beta_timecourse(ms, :, :)); % [time x trials]
                    stats.mean_beta(ms, :) = mean(beta_trace, 2);
                    stats.sem_beta(ms, :) = std(beta_trace, 0, 2) / sqrt(obj.n_trials);
                end
                
                % Baseline statistics
                baseline_data = beta_trace(baseline_idx, :);
                stats.baseline_mean(ms) = mean(baseline_data(:));
                stats.baseline_std(ms) = std(baseline_data(:));
                
                % Z-score normalization
                stats.beta_zscore(ms, :) = (stats.mean_beta(ms, :) - stats.baseline_mean(ms)) / ...
                                          stats.baseline_std(ms);
            end
            
            % Peak analysis
            stats.peak_latency = zeros(obj.n_microstates, 1);
            stats.peak_amplitude = zeros(obj.n_microstates, 1);
            
            for ms = 1:obj.n_microstates
                [max_val, max_idx] = max(abs(stats.mean_beta(ms, :)));
                stats.peak_amplitude(ms) = stats.mean_beta(ms, max_idx);
                stats.peak_latency(ms) = times(max_idx);
            end
            
            obj.display_statistics(stats, times);
        end
        
        function display_statistics(obj, stats, times)
            % Display statistical summary
            
            fprintf('\n=== Microstate Coefficient Statistics ===\n\n');
            
            fprintf('%-12s %12s %12s %12s %12s\n', ...
                    'Microstate', 'Baseline', 'Peak Amp', 'Peak Time', 'Mean |β|');
            fprintf('%s\n', repmat('-', 1, 60));
            
            for ms = 1:obj.n_microstates
                mean_abs_beta = mean(abs(stats.mean_beta(ms, :)));
                fprintf('%-12s %12.4f %12.4f %12.1f %12.4f\n', ...
                        sprintf('MS_%d', ms), ...
                        stats.baseline_mean(ms), ...
                        stats.peak_amplitude(ms), ...
                        stats.peak_latency(ms), ...
                        mean_abs_beta);
            end
            fprintf('\n');
        end
        
        function plot_coefficient_timecourses(obj, times, plot_trials)
            % Visualize coefficient time courses
            % times: time vector (ms)
            % plot_trials: boolean, plot individual trials
            
            if nargin < 2
                times = 1:obj.n_timepoints;
            end
            if nargin < 3
                plot_trials = false;
            end
            
            figure('Position', [100, 100, 1200, 800]);
            
            % Define colors for each microstate
            colors = lines(obj.n_microstates);
            
            % Plot coefficient time courses
            subplot(2, 1, 1);
            hold on;
            
            for ms = 1:obj.n_microstates
                if obj.n_trials == 1 || ~plot_trials
                    % Plot mean across trials
                    if obj.n_trials == 1
                        beta_mean = obj.beta_timecourse(ms, :);
                        beta_sem = zeros(size(beta_mean));
                    else
                        beta_data = squeeze(obj.beta_timecourse(ms, :, :));
                        beta_mean = mean(beta_data, 2);
                        beta_sem = std(beta_data, 0, 2) / sqrt(obj.n_trials);
                    end
                    
                    % Plot with shaded error
                    plot(times, beta_mean, 'Color', colors(ms, :), ...
                         'LineWidth', 2, 'DisplayName', sprintf('MS %d', ms));
                    
                    if obj.n_trials > 1
                        fill([times, fliplr(times)], ...
                             [beta_mean' + beta_sem', fliplr(beta_mean' - beta_sem')], ...
                             colors(ms, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                             'HandleVisibility', 'off');
                    end
                else
                    % Plot individual trials
                    for trial = 1:min(obj.n_trials, 10) % Limit to 10 trials
                        beta_trace = squeeze(obj.beta_timecourse(ms, :, trial));
                        plot(times, beta_trace, 'Color', [colors(ms, :), 0.3], ...
                             'LineWidth', 0.5);
                    end
                end
            end
            
            xlabel('Time (ms)');
            ylabel('Coefficient β');
            title('Microstate Coefficient Time Courses');
            legend('Location', 'best');
            grid on;
            hold off;
            
            % Plot R² time course
            subplot(2, 1, 2);
            if obj.n_trials == 1
                R2_mean = obj.R2_timecourse(:, 1);
                R2_sem = zeros(size(R2_mean));
            else
                R2_mean = mean(obj.R2_timecourse, 2);
                R2_sem = std(obj.R2_timecourse, 0, 2) / sqrt(obj.n_trials);
            end
            
            plot(times, R2_mean, 'k', 'LineWidth', 2);
            hold on;
            if obj.n_trials > 1
                fill([times, fliplr(times)], ...
                     [R2_mean' + R2_sem', fliplr(R2_mean' - R2_sem')], ...
                     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
            xlabel('Time (ms)');
            ylabel('R²');
            title('Model Fit Over Time');
            grid on;
            ylim([0, 1]);
            hold off;
        end
        
        function plot_topographies(obj, time_points, times)
            % Plot observed vs. fitted topographies at specific time points
            % time_points: indices or time values to plot
            % times: time vector (for conversion)
            
            if nargin < 3
                times = 1:obj.n_timepoints;
            end
            
            % Convert time values to indices if needed
            if max(time_points) > obj.n_timepoints
                [~, time_idx] = arrayfun(@(t) min(abs(times - t)), time_points);
            else
                time_idx = time_points;
            end
            
            n_plots = length(time_idx);
            
            % Average across trials if multiple
            if obj.n_trials == 1
                data_avg = obj.EEG;
                fitted_avg = obj.fitted_data;
                residual_avg = obj.residuals;
            else
                data_avg = mean(obj.EEG, 3);
                fitted_avg = mean(obj.fitted_data, 3);
                residual_avg = mean(obj.residuals, 3);
            end
            
            figure('Position', [100, 100, 300*n_plots, 900]);
            
            for i = 1:n_plots
                t = time_idx(i);
                
                % Observed data
                subplot(3, n_plots, i);
                if ~isempty(obj.chanlocs)
                    topoplot(data_avg(:, t), obj.chanlocs, ...
                            'maplimits', 'maxmin', 'electrodes', 'off');
                else
                    imagesc(reshape(data_avg(:, t), sqrt(obj.n_channels), []));
                    axis square;
                end
                title(sprintf('Observed\nt = %.0f ms', times(t)));
                colorbar;
                
                % Fitted data
                subplot(3, n_plots, n_plots + i);
                if ~isempty(obj.chanlocs)
                    topoplot(fitted_avg(:, t), obj.chanlocs, ...
                            'maplimits', 'maxmin', 'electrodes', 'off');
                else
                    imagesc(reshape(fitted_avg(:, t), sqrt(obj.n_channels), []));
                    axis square;
                end
                title('Fitted');
                colorbar;
                
                % Residuals
                subplot(3, n_plots, 2*n_plots + i);
                if ~isempty(obj.chanlocs)
                    topoplot(residual_avg(:, t), obj.chanlocs, ...
                            'maplimits', 'maxmin', 'electrodes', 'off');
                else
                    imagesc(reshape(residual_avg(:, t), sqrt(obj.n_channels), []));
                    axis square;
                end
                title('Residual');
                colorbar;
            end
        end
        
        function export_results(obj, filename)
            % Export results to file
            
            results.beta_timecourse = obj.beta_timecourse;
            results.fitted_data = obj.fitted_data;
            results.residuals = obj.residuals;
            results.R2_timecourse = obj.R2_timecourse;
            results.GEV = obj.GEV;
            results.microstate_maps = obj.microstate_maps;
            results.n_microstates = obj.n_microstates;
            results.n_channels = obj.n_channels;
            results.n_timepoints = obj.n_timepoints;
            results.n_trials = obj.n_trials;
            
            save(filename, 'results');
            fprintf('Results saved to: %s\n', filename);
        end
    end
end

%% Helper function: Simple LASSO solver
function beta = lasso_solver(A, y, lambda)
    % Simple coordinate descent for LASSO
    % β = argmin ||y - Aβ||² + λ||β||₁
    
    [~, p] = size(A);
    beta = zeros(p, 1);
    
    max_iter = 100;
    tol = 1e-4;
    
    for iter = 1:max_iter
        beta_old = beta;
        
        for j = 1:p
            % Compute partial residual
            r = y - A * beta + A(:, j) * beta(j);
            
            % Soft thresholding
            z = A(:, j)' * r;
            a = A(:, j)' * A(:, j);
            
            if z > lambda
                beta(j) = (z - lambda) / a;
            elseif z < -lambda
                beta(j) = (z + lambda) / a;
            else
                beta(j) = 0;
            end
        end
        
        % Check convergence
        if norm(beta - beta_old) < tol
            break;
        end
    end
end

