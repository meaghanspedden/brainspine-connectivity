function [peak_rho, peak_latency] = select_peak(t, lat_min, lat_max, zero_excl_ms)

    % Restrict latency window
    valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
    rho_vals = t(valid_idx,3);
    lag_vals = t(valid_idx,1);

    if nargin < 4 || isempty(zero_excl_ms)
        zero_excl_ms = 2; % default exclusion
    end

    % --- Find local peaks in |rho| ---
    [~, locs] = findpeaks(abs(rho_vals));

    % If no local peaks found, fall back to all points
    if isempty(locs)
        locs = (1:numel(rho_vals))';
    end

    % Exclude near-zero region
    locs = locs(abs(lag_vals(locs)) > zero_excl_ms);

    % If everything excluded, fall back to global max
    if isempty(locs)
        [~, idx_max] = max(abs(rho_vals));
        peak_rho = rho_vals(idx_max);
        peak_latency = lag_vals(idx_max);

        fprintf('All peaks within ±%g ms. Returning global max.\n', zero_excl_ms);
        return
    end

    % Sort remaining peaks by descending |rho|
    [~, ord] = sort(abs(rho_vals(locs)), 'descend');
    locs_sorted = locs(ord);

    % Select strongest as output
    idx_max = locs_sorted(1);
    peak_rho = rho_vals(idx_max);
    peak_latency = lag_vals(idx_max);

    % --- Print top 5 peaks ---
    nPrint = min(5, numel(locs_sorted));
    fprintf('\nTop %d peaks outside ±%g ms:\n', nPrint, zero_excl_ms);
    fprintf('  Rank   Lag (ms)     rho        |rho|\n');

    for k = 1:nPrint
        ii = locs_sorted(k);
        fprintf('   %d     %7.2f    %+8.4g    %8.4g\n', ...
            k, lag_vals(ii), rho_vals(ii), abs(rho_vals(ii)));
    end
end