
%% Cross-validated ridge regression for optimal regularization

function glm = fit_glm_with_cv(EEG, microstate_maps, lambda_range)
    % Fit GLM with cross-validated ridge regularization
    % lambda_range: vector of lambda values to test
    
    if nargin < 3
        lambda_range = logspace(-3, 2, 20);
    end
    
    n_folds = 5;
    [n_channels, n_timepoints, n_trials] = size(EEG);
    
    % Partition trials
    cv_indices = crossvalind('Kfold', n_trials, n_folds);
    
    % Initialize error storage
    cv_error = zeros(length(lambda_range), n_timepoints);
    
    fprintf('Cross-validation for regularization parameter...\n');
    
    for lambda_idx = 1:length(lambda_range)
        lambda = lambda_range(lambda_idx);
        fold_errors = zeros(n_folds, n_timepoints);
        
        for fold = 1:n_folds
            % Split data
            test_idx = (cv_indices == fold);
            train_idx = ~test_idx;
            
            EEG_train = EEG(:, :, train_idx);
            EEG_test = EEG(:, :, test_idx);
            
            % Fit on training data
            glm_train = MicrostateGLM(EEG_train, microstate_maps);
            glm_train.fit_glm('ridge', lambda);
            
            % Predict on test data
            for t = 1:n_timepoints
                test_data = mean(EEG_test(:, t, :), 3);
                
                % Use training coefficients
                if size(glm_train.beta_timecourse, 3) == 1
                    beta_t = glm_train.beta_timecourse(:, t);
                else
                    beta_t = mean(glm_train.beta_timecourse(:, t, :), 3);
                end
                
                predicted = microstate_maps * beta_t;
                fold_errors(fold, t) = norm(test_data - predicted)^2;
            end
        end
        
        cv_error(lambda_idx, :) = mean(fold_errors, 1);
        
        if mod(lambda_idx, 5) == 0
            fprintf('  Tested %d/%d lambda values\n', lambda_idx, length(lambda_range));
        end
    end
    
    % Find optimal lambda (minimize average error across time)
    [~, best_idx] = min(mean(cv_error, 2));
    lambda_opt = lambda_range(best_idx);
    
    fprintf('Optimal lambda: %.4f\n', lambda_opt);
    
    % Fit final model with optimal lambda
    glm = MicrostateGLM(EEG, microstate_maps);
    glm.fit_glm('ridge', lambda_opt);
    
    % Visualize CV results
    figure;
    semilogx(lambda_range, mean(cv_error, 2), 'o-', 'LineWidth', 2);
    hold on;
    plot(lambda_opt, mean(cv_error(best_idx, :)), 'r*', 'MarkerSize', 15);
    xlabel('Lambda');
    ylabel('Cross-Validation Error');
    title('Ridge Regression: Lambda Selection');
    grid on;
    legend('CV Error', 'Optimal Lambda');
end