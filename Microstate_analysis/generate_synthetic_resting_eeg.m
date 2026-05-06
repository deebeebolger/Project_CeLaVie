function EEG = generate_synthetic_resting_eeg()
    % Generate realistic synthetic resting-state EEG
    
    % Initialize
    EEG = eeg_emptyset();
    EEG.srate = 500; % 500 Hz
    duration = 120; % 120 seconds
    EEG.pnts = duration * EEG.srate;
    EEG.times = (0:EEG.pnts-1) / EEG.srate;
    
    % Load channel locations
    EEG.chanlocs = readlocs('standard_1005.elc');
    channel_subset = 1:64;
    EEG.chanlocs = EEG.chanlocs(channel_subset);
    EEG.nbchan = length(channel_subset);
    
    % Create 4 microstate templates
    n_ms = 4;
    ms_maps = create_canonical_microstates(EEG.chanlocs, n_ms);
    
    fprintf('Generating synthetic resting-state EEG...\n');
    fprintf('  Duration: %d seconds\n', duration);
    fprintf('  Sampling rate: %d Hz\n', EEG.srate);
    fprintf('  Channels: %d\n', EEG.nbchan);
    
    % Generate microstate time series with realistic dynamics
    % Use semi-Markov model for realistic transitions
    
    % Mean durations (ms)
    mean_durations = [60, 80, 90, 70];
    
    % Transition probabilities
    T = [0.0, 0.4, 0.3, 0.3;
         0.4, 0.0, 0.3, 0.3;
         0.3, 0.3, 0.0, 0.4;
         0.3, 0.4, 0.3, 0.0];
    
    % Normalize
    T = T ./ sum(T, 2);
    
    % Generate sequence
    current_ms = 1;
    t = 1;
    labels = zeros(1, EEG.pnts);
    
    while t <= EEG.pnts
        % Sample duration from gamma distribution
        duration_ms = mean_durations(current_ms);
        duration_samples = round(gamrnd(4, duration_ms/4/1000 * EEG.srate, 1, 1));
        duration_samples = max(10, min(duration_samples, EEG.pnts - t + 1));
        
        % Assign label
        labels(t:t+duration_samples-1) = current_ms;
        
        % Transition
        next_ms = randsample(n_ms, 1, true, T(current_ms, :));
        current_ms = next_ms;
        t = t + duration_samples;
    end
    
    % Generate EEG data from microstate sequence
    EEG.data = zeros(EEG.nbchan, EEG.pnts);
    
    % Add smooth transitions
    transition_samples = round(0.01 * EEG.srate); % 10 ms transition
    
    for t = 1:EEG.pnts
        current_label = labels(t);
        
        % Check for transition
        if t > 1 && labels(t) ~= labels(t-1)
            % Smooth transition
            alpha = min(1, (t - find(labels == current_label, 1, 'first')) / transition_samples);
        else
            alpha = 1;
        end
        
        % Generate topography
        topo = ms_maps(:, current_label) * alpha;
        
        % Add previous map influence during transition
        if alpha < 1 && t > 1
            prev_label = labels(t-1);
            topo = topo + ms_maps(:, prev_label) * (1 - alpha);
        end
        
        % Add oscillations (alpha rhythm)
        alpha_freq = 10; % Hz
        alpha_phase = 2 * pi * alpha_freq * EEG.times(t);
        alpha_amp = 2 + randn();
        
        % Add 1/f noise
        noise = generate_1f_noise(EEG.nbchan);
        
        EEG.data(:, t) = topo * (1 + alpha_amp * sin(alpha_phase)) + 0.5 * noise;
    end
    
    EEG = eeg_checkset(EEG);
    fprintf('Synthetic resting-state EEG generated successfully.\n');
end

function noise = generate_1f_noise(n_channels)
    % Generate 1/f (pink) noise
    
    white = randn(n_channels, 1);
    
    % Simple approximation of 1/f noise
    noise = white;
    for i = 2:n_channels
        noise(i) = 0.7 * noise(i-1) + 0.3 * white(i);
    end
end

function maps = extract_microstates_kmeans(EEG, n_microstates, gfp_peaks_only)
    % Extract microstate maps using modified K-means clustering
    
    if nargin < 3
        gfp_peaks_only = true; % Only use GFP peaks (Pascual-Marqui method)
    end
    
    fprintf('Extracting %d microstate maps using K-means...\n', n_microstates);
    
    % Prepare data
    if ndims(EEG.data) == 3
        data = reshape(EEG.data, size(EEG.data, 1), []);
    else
        data = EEG.data;
    end
    
    [n_channels, n_samples] = size(data);
    
    % Normalize each time point
    data_norm = data ./ vecnorm(data, 2, 1);
    
    % Use only GFP peaks if requested
    if gfp_peaks_only
        GFP = std(data, 0, 1);
        
        % Find local maxima
        min_distance = round(0.02 * EEG.srate); % Minimum 20 ms between peaks
        [~, peak_locs] = findpeaks(GFP, 'MinPeakDistance', min_distance);
        
        data_for_clustering = data_norm(:, peak_locs)';
        fprintf('  Using %d GFP peaks for clustering\n', length(peak_locs));
    else
        data_for_clustering = data_norm';
    end
    
    % K-means clustering with multiple restarts
    n_restarts = 20;
    best_maps = [];
    best_gev = 0;
    
    for restart = 1:n_restarts
        [idx, C] = kmeans(data_for_clustering, n_microstates, ...
                         'Distance', 'cosine', ...
                         'Replicates', 1, ...
                         'MaxIter', 1000);
        
        % Normalize cluster centers
        maps = C';
        for i = 1:n_microstates
            maps(:, i) = maps(:, i) / norm(maps(:, i));
        end
        
        % Compute GEV
        gev = compute_gev_from_maps(data_norm, maps);
        
        if gev > best_gev
            best_gev = gev;
            best_maps = maps;
        end
    end
    
    maps = best_maps;
    fprintf('  Best GEV: %.2f%%\n', best_gev * 100);
end

function gev = compute_gev_from_maps(data_norm, maps)
    % Compute Global Explained Variance
    
    n_samples = size(data_norm, 2);
    
    % Assign each sample to best matching map
    correlations = maps' * data_norm;
    [max_corr, ~] = max(abs(correlations), [], 1);
    
    gev = sum(max_corr.^2) / n_samples;
end