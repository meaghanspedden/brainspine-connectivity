function [peak_rho, peak_latency] = select_peak(t, lat_min, lat_max)
    % Restrict latency window
    valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
    rho_vals = t(valid_idx,3);
    lag_vals = t(valid_idx,1);

    % Find the global peak in absolute value
    [~, idx_max] = max(abs(rho_vals));
    peak_rho = rho_vals(idx_max);
    peak_latency = lag_vals(idx_max);

    % If the max peak is exactly zero-lag, consider secondary peak
    if peak_latency == 0
        [~, peak_locs] = findpeaks(abs(rho_vals));
        % remove zero-lag
        peak_locs(lag_vals(peak_locs) == 0) = [];
        if ~isempty(peak_locs)
            % pick largest absolute rho among remaining peaks
            [~, idx_sec] = max(abs(rho_vals(peak_locs)));
            peak_rho = rho_vals(peak_locs(idx_sec));
            peak_latency = lag_vals(peak_locs(idx_sec));
        end
    end
end