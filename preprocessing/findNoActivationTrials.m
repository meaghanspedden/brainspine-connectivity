function no_activation_indices = findNoActivationTrials(EMG_data, threshold)
    % EMG_data: A matrix where each column represents a trial and each row represents time or samples
    % threshold: The variance threshold below which trials will be considered as 'no activation'

    % Calculate the variance for each trial
    trial_variance = var(EMG_data, 0, 1);  % 0 for population variance, 1 for calculating variance along rows (each trial)

    % Find trials where variance is below the threshold
    no_activation_indices = find(trial_variance < threshold);
end
