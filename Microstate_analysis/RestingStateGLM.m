
%% Spatial GLM for Resting-State EEG Analysis
% Models continuous resting-state EEG as a linear combination of spatial templates
% Extracts microstate dynamics, transitions, and network properties

classdef RestingStateGLM < handle
    properties
        EEG_data            % [n_channels x n_timepoints] continuous resting EEG
        microstate_maps     % [n_channels x n_microstates] spatial templates
        chanlocs            % Channel locations
        srate               % Sampling rate (Hz)
        
        % GLM results
        beta_timecourse     % [n_microstates x n_timepoints] coefficient time series
        fitted_data         % Reconstructed EEG
        residuals           % Residual signal
        
        % Microstate metrics
        MS_labels           % [1 x n_timepoints] assigned microstate at each sample
        MS_parameters       % Structure with microstate parameters
        
        % Spatial statistics
        W                   % Spatial weights matrix
        spatial_stats       % Spatial autocorrelation results
        
        % Model quality
        GEV                 % Global Explained Variance
        R2_timecourse       % Time-varying R²
    end
    
    methods
        function obj = RestingStateGLM(EEG, microstate_maps, chanlocs, srate)
            % Constructor
            % EEG: [n_channels x n_timepoints] or EEGLAB structure
            % microstate_maps: [n_channels x n_microstates]
            % chanlocs: channel locations structure
            % srate: sampling rate
            
            if isstruct(EEG)
                obj.EEG_data = double(EEG.data(:, :));
                obj.chanlocs = EEG.chanlocs;
                obj.srate = EEG.srate;
            else
                obj.EEG_data = double(EEG);
                obj.chanlocs = chanlocs;
                obj.srate = srate;
            end
            
            % Remove any 3D dimension (trials)
            if ndims(obj.EEG_data) == 3
                obj.EEG_data = reshape(obj.EEG_data, size(obj.EEG_data, 1), []);
            end
            
            obj.microstate_maps = double(microstate_maps);
            
            % Normalize microstate maps
            for i = 1:size(microstate_maps, 2)
                obj.microstate_maps(:, i) = obj.microstate_maps(:, i) / ...
                    norm(obj.microstate_maps(:, i));
            end
            
            fprintf('Initialized Resting-State GLM:\n');
            fprintf('  Channels: %d\n', size(obj.EEG_data, 1));
            fprintf('  Samples: %d (%.1f seconds)\n', size(obj.EEG_data, 2), ...
                    size(obj.EEG_data, 2) / obj.srate);
            fprintf('  Microstates: %d\n', size(obj.microstate_maps, 2));
        end
        
        function fit_continuous_glm(obj, method, params)
            % Fit GLM to continuous resting-state data
            % method: 'ols', 'ridge', 'robust', 'nnls', 'tapered_sliding'
            % params: structure with method-specific parameters
            
            if nargin < 2, method = 'ols'; end
            if nargin < 3, params = struct(); end
            
            [n_channels, n_timepoints] = size(obj.EEG_data);
            n_microstates = size(obj.microstate_maps, 2);
            
            fprintf('\nFitting continuous GLM using %s method...\n', upper(method));
            
            % Initialize
            obj.beta_timecourse = zeros(n_microstates, n_timepoints);
            obj.fitted_data = zeros(n_channels, n_timepoints);
            obj.R2_timecourse = zeros(1, n_timepoints);
            
            A = obj.microstate_maps; % Design matrix
            
            switch lower(method)
                case 'ols'
                    % Standard OLS at each time point
                    obj.fit_ols_pointwise(A);
                    
                case 'tapered_sliding'
                    % Tapered sliding window for temporal smoothness
                    if ~isfield(params, 'window_size')
                        params.window_size = round(0.1 * obj.srate); % 100ms default
                    end
                    obj.fit_tapered_sliding_window(A, params.window_size);
                    
                case 'ridge'
                    % Ridge regression with spatial regularization
                    if ~isfield(params, 'lambda')
                        params.lambda = 0.01;
                    end
                    obj.fit_ridge_pointwise(A, params.lambda);
                    
                case 'robust'
                    % Robust regression (resistant to outliers)
                    obj.fit_robust_pointwise(A);
                    
                case 'nnls'
                    % Non-negative least squares
                    obj.fit_nnls_pointwise(A);
                    
                case 'aahc'
                    % Atomize and Agglomerate Hierarchical Clustering
                    % (Modified Pascual-Marqui method)
                    if ~isfield(params, 'smoothing')
                        params.smoothing = 5; % samples
                    end
                    obj.fit_aahc_method(A, params.smoothing);
            end
            
            % Compute residuals
            obj.residuals = obj.EEG_data - obj.fitted_data;
            
            % Global metrics
            obj.compute_global_metrics();
            
            fprintf('  GEV: %.2f%%\n', obj.GEV * 100);
            fprintf('  Mean R²: %.4f\n', mean(obj.R2_timecourse));
        end
        
        function fit_ols_pointwise(obj, A)
            % Standard OLS at each time point independently
            
            [n_channels, n_timepoints] = size(obj.EEG_data);
            ATA_inv = (A' * A) \ eye(size(A, 2));
            
            for t = 1:n_timepoints
                y = obj.EEG_data(:, t);
                
                % β = (A'A)^(-1) A'y
                beta = ATA_inv * (A' * y);
                obj.beta_timecourse(:, t) = beta;
                
                % Fitted value
                y_fit = A * beta;
                obj.fitted_data(:, t) = y_fit;
                
                % R²
                SS_res = sum((y - y_fit).^2);
                SS_tot = sum((y - mean(y)).^2);
                obj.R2_timecourse(t) = max(0, 1 - SS_res / SS_tot);
                
                if mod(t, 10000) == 0
                    fprintf('    Processed %d/%d samples\n', t, n_timepoints);
                end
            end
        end
        
        function fit_tapered_sliding_window(obj, A, window_size)
            % Fit GLM using tapered sliding window for temporal smoothness
            % Reduces noise while preserving rapid transitions
            
            [n_channels, n_timepoints] = size(obj.EEG_data);
            half_win = floor(window_size / 2);
            
            % Create Hanning window for tapering
            hann_win = hanning(window_size);
            hann_win = hann_win / sum(hann_win);
            
            fprintf('    Using tapered window: %d samples (%.1f ms)\n', ...
                    window_size, window_size * 1000 / obj.srate);
            
            for t = 1:n_timepoints
                % Define window bounds
                t_start = max(1, t - half_win);
                t_end = min(n_timepoints, t + half_win);
                
                % Extract windowed data
                y_window = obj.EEG_data(:, t_start:t_end);
                
                % Create appropriate window weights
                win_idx = (t_start:t_end) - t + half_win + 1;
                win_idx = max(1, min(window_size, win_idx));
                weights = hann_win(win_idx);
                
                % Weighted regression
                y_weighted = y_window .* weights';
                y_center = sum(y_weighted, 2);
                
                % Solve
                beta = (A' * A) \ (A' * y_center);
                obj.beta_timecourse(:, t) = beta;
                
                % Fitted value at center point
                y_fit = A * beta;
                obj.fitted_data(:, t) = y_fit;
                
                % R²
                y_actual = obj.EEG_data(:, t);
                SS_res = sum((y_actual - y_fit).^2);
                SS_tot = sum((y_actual - mean(y_actual)).^2);
                obj.R2_timecourse(t) = max(0, 1 - SS_res / SS_tot);
                
                if mod(t, 10000) == 0
                    fprintf('    Processed %d/%d samples\n', t, n_timepoints);
                end
            end
        end
        
        function fit_ridge_pointwise(obj, A, lambda)
            % Ridge regression at each time point
            
            [n_channels, n_timepoints] = size(obj.EEG_data);
            n_ms = size(A, 2);
            
            % Precompute (A'A + λI)^(-1) A'
            ridge_matrix = (A' * A + lambda * eye(n_ms)) \ A';
            
            for t = 1:n_timepoints
                y = obj.EEG_data(:, t);
                beta = ridge_matrix * y;
                
                obj.beta_timecourse(:, t) = beta;
                obj.fitted_data(:, t) = A * beta;
                
                SS_res = sum((y - obj.fitted_data(:, t)).^2);
                SS_tot = sum((y - mean(y)).^2);
                obj.R2_timecourse(t) = max(0, 1 - SS_res / SS_tot);
                
                if mod(t, 10000) == 0
                    fprintf('    Processed %d/%d samples\n', t, n_timepoints);
                end
            end
        end
        
        function fit_robust_pointwise(obj, A)
            % Robust regression at each time point
            
            [~, n_timepoints] = size(obj.EEG_data);
            
            for t = 1:n_timepoints
                y = obj.EEG_data(:, t);
                
                % Robust regression using iteratively reweighted LS
                beta = robustfit(A, y, 'huber', [], 'off');
                obj.beta_timecourse(:, t) = beta;
                obj.fitted_data(:, t) = A * beta;
                
                SS_res = sum((y - obj.fitted_data(:, t)).^2);
                SS_tot = sum((y - mean(y)).^2);
                obj.R2_timecourse(t) = max(0, 1 - SS_res / SS_tot);
                
                if mod(t, 5000) == 0
                    fprintf('    Processed %d/%d samples\n', t, n_timepoints);
                end
            end
        end
        
        function fit_nnls_pointwise(obj, A)
            % Non-negative least squares at each time point
            
            [~, n_timepoints] = size(obj.EEG_data);
            
            for t = 1:n_timepoints
                y = obj.EEG_data(:, t);
                
                % NNLS: min ||Aβ - y||² subject to β ≥ 0
                beta = lsqnonneg(A, y);
                obj.beta_timecourse(:, t) = beta;
                obj.fitted_data(:, t) = A * beta;
                
                SS_res = sum((y - obj.fitted_data(:, t)).^2);
                SS_tot = sum((y - mean(y)).^2);
                obj.R2_timecourse(t) = max(0, 1 - SS_res / SS_tot);
                
                if mod(t, 10000) == 0
                    fprintf('    Processed %d/%d samples\n', t, n_timepoints);
                end
            end
        end
        
        function fit_aahc_method(obj, A, smoothing)
            % Modified AAHC (Atomize and Agglomerate Hierarchical Clustering)
            % Assigns discrete microstate labels with smooth transitions