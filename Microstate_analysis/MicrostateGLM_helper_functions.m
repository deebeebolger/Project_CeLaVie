%% Create synthetic ERP data

function EEG = create_synthetic_erp_data()
    % Create realistic synthetic ERP data with microstate structure
    
    % Initialize EEGLAB structure
    EEG = eeg_emptyset();
    
    % Load standard channel locations
    EEG.chanlocs = readlocs('standard_1005.elc');
    channel_subset = 1:64;
    EEG.chanlocs = EEG.chanlocs(channel_subset);
    EEG.nbchan = length(channel_subset);
    
    % Time parameters
    EEG.srate = 500; % 500 Hz
    EEG.times = -200:2:800; % -200 to 800 ms
    EEG.pnts = length(EEG.times);
    EEG.trials = 50;
    
    % Create 4 microstate templates
    n_ms = 4;
    ms_maps = create_canonical_microstates(EEG.chanlocs, n_ms);
    
    % Generate time courses for each microstate with ERP structure
    t = EEG.times / 1000; % Convert to seconds
    
    % Define microstate activations over time
    % MS1: Early sensory response (50-150ms)
    ms1_time = exp(-((t - 0.1).^2) / (2 * 0.02^2));
    
    % MS2: Mid-latency component (150-250ms)
    ms2_time = exp(-((t - 0.2).^2) / (2 * 0.03^2));
    
    % MS3: Late positive component (250-400ms)
    ms3_time = exp(-((t - 0.35).^2) / (2 * 0.05^2));
    
    % MS4: Sustained activity (400-600ms)
    ms4_time = (t > 0.4 & t < 0.6) .* (1 - abs(t - 0.5) / 0.1);
    
    % Combine into time course matrix
    beta_true = [ms1_time; ms2_time; ms3_time; ms4_time];
    
    % Generate EEG data
    EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
    
    for trial = 1:EEG.trials
        % Add trial-to-trial variability
        trial_beta = beta_true .* (1 + 0.3 * randn(n_ms, 1));
        
        % Synthesize EEG: data = Σ β_i(t) * map_i + noise
        for t_idx = 1:EEG.pnts
            eeg_t = ms_maps * trial_beta(:, t_idx);
            noise = 0.5 * randn(EEG.nbchan, 1);
            EEG.data(:, t_idx, trial) = eeg_t + noise;
        end
    end
    
    EEG = eeg_checkset(EEG);
    
    fprintf('Created synthetic ERP data:\n');
    fprintf('  %d channels, %d time points, %d trials\n', ...
            EEG.nbchan, EEG.pnts, EEG.trials);
end

%% Create canonical microstate maps
function maps = create_canonical_microstates(chanlocs, n_microstates)
    % Create canonical microstate topographies
    % Based on Koenig et al. (2002) patterns
    
    n_channels = length(chanlocs);
    maps = zeros(n_channels, n_microstates);
    
    % Extract channel positions
    theta = [chanlocs.theta]' * pi / 180; % Convert to radians
    radius = [chanlocs.radius]';
    
    % Normalize radius
    radius = radius / max(radius);
    
    % Convert to Cartesian
    x = radius .* cos(theta);
    y = radius .* sin(theta);
    
    % Define microstate patterns
    for ms = 1:n_microstates
        switch ms
            case 1
                % Left-right pattern
                maps(:, ms) = x;
                
            case 2
                % Anterior-posterior pattern
                maps(:, ms) = y;
                
            case 3
                % Diagonal pattern
                maps(:, ms) = x + y;
                
            case 4
                % Radial pattern
                maps(:, ms) = radius .* cos(2 * theta);
                
            otherwise
                % Random pattern
                maps(:, ms) = randn(n_channels, 1);
        end
        
        % Normalize
        maps(:, ms) = maps(:, ms) / norm(maps(:, ms));
    end
    
    % Add some spatial smoothness
    for ms = 1:n_microstates
        % Simple spatial smoothing using distance
        coords = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'];
        D = pdist2(coords, coords);
        W = exp(-D.^2 / (2 * 0.1^2));
        W = W ./ sum(W, 2);
        
        maps(:, ms) = W * maps(:, ms);
        maps(:, ms) = maps(:, ms) / norm(maps(:, ms));
    end
end

%% Cluster-based permutation test
function results = cluster_permutation_test(glm, baseline_window, times, n_perm)
    % Cluster-based permutation test for coefficient time courses
    % Tests: H0: β(t) = baseline for each microstate
    
    if nargin < 4, n_perm = 1000; end
    
    fprintf('Running %d permutations...\n', n_perm);
    
    % Find baseline indices
    baseline_idx = times >= baseline_window(1) & times <= baseline_window(2);
    
    n_ms = glm.n_microstates;
    n_time = glm.n_timepoints;
    
    % Compute observed t-statistics
    if glm.n_trials == 1
        error('Permutation test requires multiple trials');
    end
    
    observed_t = zeros(n_ms, n_time);
    cluster_stats = cell(n_ms, 1);
    
    for ms = 1:n_ms
        beta_data = squeeze(glm.beta_timecourse(ms, :, :))'; % [trials x time]
        baseline_mean = mean(beta_data(:, baseline_idx), 2);
        
        % Subtract baseline from each trial
        beta_corrected = beta_data - baseline_mean;
        
        % T-test at each time point
        [~, ~, ~, stats] = ttest(beta_corrected);
        observed_t(ms, :) = stats.tstat;
        
        % Find clusters (contiguous significant time points)
        sig_mask = abs(observed_t(ms, :)) > tinv(0.975, glm.n_trials - 1);
        cluster_stats{ms} = find_clusters(sig_mask, observed_t(ms, :));
    end
    
    % Permutation distribution
    max_cluster_mass = zeros(n_perm, n_ms);
    
    for perm = 1:n_perm
        for ms = 1:n_ms
            beta_data = squeeze(glm.beta_timecourse(ms, :, :))';
            baseline_mean = mean(beta_data(:, baseline_idx), 2);
            beta_corrected = beta_data - baseline_mean;
            
            % Random sign flip
            signs = 2 * (rand(glm.n_trials, 1) > 0.5) - 1;
            beta_perm = beta_corrected .* signs;
            
            % Compute permuted t-statistics
            [~, ~, ~, stats] = ttest(beta_perm);
            t_perm = stats.tstat;
            
            % Find clusters
            sig_mask = abs(t_perm) > tinv(0.975, glm.n_trials - 1);
            clusters_perm = find_clusters(sig_mask, t_perm);
            
            % Store maximum cluster mass
            if ~isempty(clusters_perm)
                max_cluster_mass(perm, ms) = max([clusters_perm.mass]);
            end
        end
        
        if mod(perm, 100) == 0
            fprintf('  Permutation %d/%d\n', perm, n_perm);
        end
    end
    
    % Compute p-values for observed clusters
    results.observed_t = observed_t;
    results.cluster_stats = cluster_stats;
    results.max_cluster_mass = max_cluster_mass;
    results.alpha = 0.05;
    
    for ms = 1:n_ms
        if ~isempty(cluster_stats{ms})
            for c = 1:length(cluster_stats{ms})
                cluster_mass = cluster_stats{ms}(c).mass;
                p_value = mean(max_cluster_mass(:, ms) >= cluster_mass);
                cluster_stats{ms}(c).p_value = p_value;
                
                fprintf('MS%d Cluster %d: t_sum=%.2f, p=%.4f\n', ...
                        ms, c, cluster_mass, p_value);
            end
        end
    end
    
    results.cluster_stats = cluster_stats;
end

%% Find temporal clusters
function clusters = find_clusters(sig_mask, t_values)
    % Find contiguous clusters of significant time points
    
    clusters = [];
    
    if ~any(sig_mask)
        return;
    end
    
    % Find cluster boundaries
    diff_mask = diff([0, sig_mask, 0]);
    starts = find(diff_mask == 1);
    ends = find(diff_mask == -1) - 1;
    
    % Compute cluster statistics
    for i = 1:length(starts)
        cluster.start = starts(i);
        cluster.end = ends(i);
        cluster.length = ends(i) - starts(i) + 1;
        cluster.mass = sum(abs(t_values(starts(i):ends(i))));
        clusters = [clusters, cluster];
    end
end