%% Spatial GLM for EEG Scalp Maps
% This code implements a spatial general linear model where EEG topographic
% maps at different time points or conditions serve as spatial regressors

classdef SpatialGLM_EEG < handle
    properties
        chanlocs        % Channel locations structure
        n_channels      % Number of EEG channels
        W               % Spatial weights matrix
        design_matrix   % GLM design matrix (spatial maps)
        dependent_var   % Dependent variable to predict
        beta            % GLM coefficients
        residuals       % Model residuals
        fitted_values   % Fitted values
        statistics      % Statistical results
    end
    
    methods
        function obj = SpatialGLM_EEG(chanlocs)
            % Constructor
            obj.chanlocs = chanlocs;
            obj.n_channels = length(chanlocs);
            obj.create_spatial_weights();
        end
        
        function create_spatial_weights(obj, method, k)
            % Create spatial weights matrix based on electrode positions
            % method: 'knn' (default), 'distance', 'inverse_distance'
            % k: number of neighbors or distance threshold
            
            if nargin < 2, method = 'knn'; end
            if nargin < 3, k = 4; end
            
            % Extract 3D electrode coordinates
            coords = zeros(obj.n_channels, 3);
            for i = 1:obj.n_channels
                coords(i, :) = [obj.chanlocs(i).X, ...
                               obj.chanlocs(i).Y, ...
                               obj.chanlocs(i).Z];
            end
            
            % Handle 2D projection if Z is missing
            if all(coords(:,3) == 0)
                coords = coords(:, 1:2);
            end
            
            % Create weights matrix
            n = obj.n_channels;
            W = zeros(n, n);
            
            switch method
                case 'knn'
                    % K-nearest neighbors
                    [idx, dist] = knnsearch(coords, coords, 'K', k+1);
                    idx = idx(:, 2:end); % Remove self
                    
                    for i = 1:n
                        W(i, idx(i,:)) = 1;
                    end
                    
                case 'distance'
                    % Distance threshold
                    D = pdist2(coords, coords);
                    W = (D <= k & D > 0);
                    
                case 'inverse_distance'
                    % Inverse distance weighting
                    D = pdist2(coords, coords);
                    W = 1 ./ D;
                    W(isinf(W)) = 0;
                    % Apply distance threshold
                    if k > 0
                        W(D > k) = 0;
                    end
                    
                case 'gaussian'
                    % Gaussian kernel
                    D = pdist2(coords, coords);
                    sigma = k; % bandwidth parameter
                    W = exp(-D.^2 / (2 * sigma^2));
                    W(eye(n) == 1) = 0; % Remove self-connections
            end
            
            % Row-standardize
            row_sums = sum(W, 2);
            row_sums(row_sums == 0) = 1; % Avoid division by zero
            W = W ./ row_sums;
            
            obj.W = W;
        end
        
        function build_design_matrix(obj, scalp_maps, include_spatial_lag)
            % Build GLM design matrix from EEG scalp maps
            % scalp_maps: [n_channels x n_maps] matrix
            % Each column is a topographic map serving as a regressor
            % include_spatial_lag: boolean, add spatially lagged regressors
            
            if nargin < 3, include_spatial_lag = false; end
            
            % Normalize each scalp map
            n_maps = size(scalp_maps, 2);
            X = zeros(obj.n_channels, n_maps);
            
            for i = 1:n_maps
                map = scalp_maps(:, i);
                % Z-score normalization
                X(:, i) = (map - mean(map)) / std(map);
            end
            
            % Add constant term
            X = [ones(obj.n_channels, 1), X];
            
            % Optionally add spatially lagged regressors (WX)
            if include_spatial_lag
                WX = obj.W * scalp_maps;
                % Normalize lagged maps
                for i = 1:size(WX, 2)
                    WX(:, i) = (WX(:, i) - mean(WX(:, i))) / std(WX(:, i));
                end
                X = [X, WX];
            end
            
            obj.design_matrix = X;
        end
        
        function fit_ols(obj, y)
            % Fit standard OLS GLM
            % y: dependent variable [n_channels x 1]
            
            obj.dependent_var = y;
            X = obj.design_matrix;
            
            % OLS estimation: β = (X'X)^(-1) X'y
            obj.beta = (X' * X) \ (X' * y);
            obj.fitted_values = X * obj.beta;
            obj.residuals = y - obj.fitted_values;
            
            % Compute statistics
            n = obj.n_channels;
            k = size(X, 2);
            
            % Residual sum of squares
            RSS = sum(obj.residuals.^2);
            TSS = sum((y - mean(y)).^2);
            
            % R-squared
            R2 = 1 - RSS / TSS;
            R2_adj = 1 - (1 - R2) * (n - 1) / (n - k);
            
            % Standard errors
            sigma2 = RSS / (n - k);
            var_beta = sigma2 * inv(X' * X);
            se_beta = sqrt(diag(var_beta));
            
            % T-statistics
            t_stats = obj.beta ./ se_beta;
            p_values = 2 * (1 - tcdf(abs(t_stats), n - k));
            
            % Store statistics
            obj.statistics.method = 'OLS';
            obj.statistics.R2 = R2;
            obj.statistics.R2_adj = R2_adj;
            obj.statistics.sigma2 = sigma2;
            obj.statistics.se_beta = se_beta;
            obj.statistics.t_stats = t_stats;
            obj.statistics.p_values = p_values;
            obj.statistics.RSS = RSS;
            obj.statistics.TSS = TSS;
            
            % Compute AIC and BIC
            loglik = -n/2 * log(2*pi) - n/2 * log(sigma2) - RSS/(2*sigma2);
            obj.statistics.loglik = loglik;
            obj.statistics.AIC = -2 * loglik + 2 * k;
            obj.statistics.BIC = -2 * loglik + log(n) * k;
        end
        
        function fit_spatial_lag(obj, y)
            % Fit Spatial Lag Model: y = ρWy + Xβ + ε
            % y: dependent variable [n_channels x 1]
            
            obj.dependent_var = y;
            X = obj.design_matrix;
            W = obj.W;
            n = obj.n_channels;
            k = size(X, 2);
            
            % Concentrated log-likelihood
            function loglik = conc_loglik(rho)
                eig_W = eig(full(W));
                log_det = sum(log(1 - rho * eig_W));
                
                y_star = y - rho * W * y;
                beta = X \ y_star;
                e = y_star - X * beta;
                
                sigma2 = (e' * e) / n;
                loglik = -0.5*n*log(2*pi) - 0.5*n*log(sigma2) + log_det;
                loglik = -loglik; % Minimize negative
            end
            
            % Optimize
            eig_W = eig(full(W));
            rho_bounds = [1/min(eig_W), 1/max(eig_W)];
            options = optimset('Display', 'off', 'TolX', 1e-6);
            
            [rho_hat, neg_loglik] = fminbnd(@conc_loglik, ...
                                            rho_bounds(1)*0.99, ...
                                            rho_bounds(2)*0.99, ...
                                            options);
            
            % Final estimates
            y_star = y - rho_hat * W * y;
            beta_hat = X \ y_star;
            obj.beta = beta_hat;
            obj.fitted_values = y - (y_star - X * beta_hat);
            obj.residuals = y - obj.fitted_values;
            
            % Statistics
            RSS = sum(obj.residuals.^2);
            TSS = sum((y - mean(y)).^2);
            sigma2_hat = (obj.residuals' * obj.residuals) / n;
            
            % Variance-covariance (simplified)
            var_beta = sigma2_hat * inv(X' * X);
            se_beta = sqrt(diag(var_beta));
            
            t_stats = beta_hat ./ se_beta;
            p_values = 2 * (1 - tcdf(abs(t_stats), n - k - 1));
            
            % Store results
            obj.statistics.method = 'Spatial Lag';
            obj.statistics.rho = rho_hat;
            obj.statistics.R2 = 1 - RSS / TSS;
            obj.statistics.sigma2 = sigma2_hat;
            obj.statistics.se_beta = se_beta;
            obj.statistics.t_stats = t_stats;
            obj.statistics.p_values = p_values;
            obj.statistics.loglik = -neg_loglik;
            obj.statistics.AIC = -2*(-neg_loglik) + 2*(k+1);
            obj.statistics.BIC = -2*(-neg_loglik) + log(n)*(k+1);
        end
        
        function [I, p_value] = test_spatial_autocorr(obj, use_residuals)
            % Moran's I test for spatial autocorrelation
            % use_residuals: test residuals (true) or dependent var (false)
            
            if nargin < 2, use_residuals = true; end
            
            if use_residuals && ~isempty(obj.residuals)
                data = obj.residuals;
            else
                data = obj.dependent_var;
            end
            
            [I, p_value] = morans_i_eeg(data, obj.W);
        end
        
        function display_results(obj)
            % Display GLM results
            
            fprintf('\n========================================\n');
            fprintf('  Spatial GLM Results (%s)\n', obj.statistics.method);
            fprintf('========================================\n\n');
            
            fprintf('Model Fit:\n');
            if isfield(obj.statistics, 'R2')
                fprintf('  R² = %.4f\n', obj.statistics.R2);
            end
            if isfield(obj.statistics, 'R2_adj')
                fprintf('  Adjusted R² = %.4f\n', obj.statistics.R2_adj);
            end
            fprintf('  AIC = %.2f\n', obj.statistics.AIC);
            fprintf('  BIC = %.2f\n', obj.statistics.BIC);
            fprintf('  Log-likelihood = %.2f\n', obj.statistics.loglik);
            
            if isfield(obj.statistics, 'rho')
                fprintf('\nSpatial Parameter:\n');
                fprintf('  ρ (rho) = %.4f\n', obj.statistics.rho);
            end
            
            fprintf('\nCoefficients:\n');
            fprintf('%-15s %10s %10s %10s %10s\n', ...
                    'Predictor', 'Beta', 'SE', 't-stat', 'p-value');
            fprintf('%-15s %10s %10s %10s %10s\n', ...
                    repmat('-', 1, 15), repmat('-', 1, 10), ...
                    repmat('-', 1, 10), repmat('-', 1, 10), ...
                    repmat('-', 1, 10));
            
            n_pred = length(obj.beta);
            pred_names = cell(n_pred, 1);
            pred_names{1} = 'Intercept';
            for i = 2:n_pred
                pred_names{i} = sprintf('Predictor_%d', i-1);
            end
            
            for i = 1:n_pred
                fprintf('%-15s %10.4f %10.4f %10.4f %10.4f\n', ...
                        pred_names{i}, obj.beta(i), ...
                        obj.statistics.se_beta(i), ...
                        obj.statistics.t_stats(i), ...
                        obj.statistics.p_values(i));
            end
            fprintf('\n');
        end
        
        function plot_results(obj)
            % Visualize GLM results on scalp topography
            
            figure('Position', [100, 100, 1400, 400]);
            
            % Plot 1: Dependent variable
            subplot(1, 4, 1);
            topoplot(obj.dependent_var, obj.chanlocs, ...
                    'maplimits', 'maxmin', 'electrodes', 'on');
            title('Dependent Variable');
            colorbar;
            
            % Plot 2: Fitted values
            subplot(1, 4, 2);
            topoplot(obj.fitted_values, obj.chanlocs, ...
                    'maplimits', 'maxmin', 'electrodes', 'on');
            title('Fitted Values');
            colorbar;
            
            % Plot 3: Residuals
            subplot(1, 4, 3);
            topoplot(obj.residuals, obj.chanlocs, ...
                    'maplimits', 'maxmin', 'electrodes', 'on');
            title('Residuals');
            colorbar;
            
            % Plot 4: Residuals vs Fitted
            subplot(1, 4, 4);
            scatter(obj.fitted_values, obj.residuals, 30, 'filled');
            xlabel('Fitted Values');
            ylabel('Residuals');
            title('Residuals vs Fitted');
            grid on;
            hold on;
            plot(xlim, [0 0], 'r--', 'LineWidth', 1.5);
        end
    end
end

%% Helper function: Moran's I for EEG
function [I, p_value] = morans_i_eeg(data, W)
    % Moran's I statistic for spatial autocorrelation
    
    n = length(data);
    data = data - mean(data);
    
    numerator = data' * W * data;
    denominator = data' * data;
    S0 = sum(W(:));
    
    I = (n / S0) * (numerator / denominator);
    
    % Expected value and variance
    E_I = -1 / (n - 1);
    
    % Variance calculation
    S1 = 0.5 * sum((W + W').^2, 'all');
    S2 = sum((sum(W, 2) + sum(W, 1)').^2);
    
    Var_I = ((n * S1 - 2 * S2) / ((n-1)*(n-2)*S0^2)) - E_I^2;
    
    % Z-score and p-value
    z = (I - E_I) / sqrt(Var_I);
    p_value = 2 * (1 - normcdf(abs(z)));
    
    fprintf('Moran''s I = %.4f (p = %.4f)\n', I, p_value);
end