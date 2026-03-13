%% get_brain_maxima_mni.m
%
% Standalone script: loads saved brain group results from Steps 1 and 4
% and reports MNI coordinates of:
%   (a) Participant 1 single-subject maximum (native -> MNI via T)
%   (b) Group prevalence maximum (native -> MNI via T)
%   (c) All tied/near-tied maxima within a relative tolerance
%
% The four maps covered:
%   Step 1:  Brain-EMG coherence        (P1 + group)
%   Step 4:  Brain-spineVE coherence    (P1 + group)
%
% EDIT USER CONFIG below, then run.

%% ---- USER CONFIG --------------------------------------------------------
save_dir       = 'C:\Users\mspedden\Documents\brainspine_save_bemv2';
T_mat_path     = 'C:\Users\mspedden\Documents\brainspine_save\T.mat';
geomfile       = 'C:\Users\mspedden\Documents\new_leadfields_and_geom\geometries_cervical_realistic.mat';
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';

% filenames — adjust suffix to match what you ran
step1_file = 'groupRes_brain_DICS_bemv2_brainEMG_brainSmooth_8mm.mat';
step4_file = 'groupRes_brain_DICS_spineVC_bemv2_functionalVE_spineSmooth_20mm_brainSmooth_8mm.mat';

% prevalence threshold (must match pipeline)
brain_prev_threshold = 0.3;

% RELATIVE tolerance: sources within this *fraction* of the global max are
% also reported as peaks. E.g. 0.01 = within 1% of max. Set to 0 for
% single peak only (exact max, which may have ties on a discrete grid).
peak_tol_frac = 0.02;

% participant index for "Participant 1" in the subjResults array
p1_idx = 1;
%% -----------------------------------------------------------------------

addpath(fieldtrip_path);
ft_defaults;

%% load T matrix
load(T_mat_path, 'T');
T_inv = inv(T);

%% =========================================================================
%  STEP 1 -- Brain-EMG coherence
%% =========================================================================
fprintf('\n=================================================================\n');
fprintf(' STEP 1 -- Brain-EMG coherence\n');
fprintf('=================================================================\n');

s1  = load(fullfile(save_dir, step1_file), 'subjResults');
sr1 = s1.subjResults;

% --- P1: use stored maxpos_mni (computed in pipeline from max coh_diff) ---
fprintf('\n  Participant 1 maximum (stored in pipeline output):\n');
mni_p1 = sr1(p1_idx).maxpos_mni;
for k = 1:size(mni_p1,1)
    fprintf('    Peak %d:  MNI [%.0f  %.0f  %.0f]\n', ...
        k, mni_p1(k,1), mni_p1(k,2), mni_p1(k,3));
end

% --- P1 cross-check: recompute from coh_diff with relative tolerance ---
pos_inside1 = sr1(p1_idx).pos(logical(sr1(p1_idx).inside), :);
report_peaks('Participant 1 -- coh_diff cross-check', ...
    sr1(p1_idx).coh_diff(:), pos_inside1, T_inv, peak_tol_frac);

% --- Group prevalence maximum ---
all_masks1   = cat(2, sr1(:).sig_mask);
prev1        = mean(all_masks1, 2);
prev1_masked = prev1;
prev1_masked(prev1_masked < brain_prev_threshold) = 0;
pos_all1     = sr1(1).pos;
inside1      = logical(sr1(1).inside);

report_peaks('Group prevalence maximum -- Brain-EMG', ...
    prev1_masked(inside1), pos_all1(inside1,:), T_inv, peak_tol_frac);

%% =========================================================================
%  STEP 4 -- Brain-spineVE coherence
%% =========================================================================
fprintf('\n=================================================================\n');
fprintf(' STEP 4 -- Brain-spineVE coherence\n');
fprintf('=================================================================\n');

s4  = load(fullfile(save_dir, step4_file), 'subjResults');
sr4 = s4.subjResults;

% --- P1: use stored maxpos_mni ---
fprintf('\n  Participant 1 maximum (stored in pipeline output):\n');
mni_p4 = sr4(p1_idx).maxpos_mni;
for k = 1:size(mni_p4,1)
    fprintf('    Peak %d:  MNI [%.0f  %.0f  %.0f]\n', ...
        k, mni_p4(k,1), mni_p4(k,2), mni_p4(k,3));
end

% --- P1 cross-check ---
pos_inside4 = sr4(p1_idx).pos(logical(sr4(p1_idx).inside), :);
report_peaks('Participant 1 -- coh_diff cross-check', ...
    sr4(p1_idx).coh_diff(:), pos_inside4, T_inv, peak_tol_frac);

% --- Group prevalence maximum ---
all_masks4   = cat(2, sr4(:).sig_mask);
prev4        = mean(all_masks4, 2);
prev4_masked = prev4;
prev4_masked(prev4_masked < brain_prev_threshold) = 0;
pos_all4     = sr4(1).pos;
inside4      = logical(sr4(1).inside);

report_peaks('Group prevalence maximum -- Brain-spineVE', ...
    prev4_masked(inside4), pos_all4(inside4,:), T_inv, peak_tol_frac);

fprintf('\nDone.\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function mni = native_to_mni(pos_mm, T_inv)
    ph  = [pos_mm, ones(size(pos_mm,1),1)]';
    mni = (T_inv * ph)';
    mni = mni(:,1:3);
end

function report_peaks(label, vals, pos_native, T_inv, tol_frac)
    vals   = vals(:);
    maxval = max(vals);

    % relative tolerance: threshold = max * (1 - tol_frac)
    thresh    = maxval * (1 - tol_frac);
    peak_mask = vals >= thresh;
    peak_pos  = pos_native(peak_mask, :);
    peak_vals = vals(peak_mask);
    peak_mni  = native_to_mni(peak_pos, T_inv);

    % sort by value descending
    [peak_vals, ord] = sort(peak_vals, 'descend');
    peak_pos = peak_pos(ord,:);
    peak_mni = peak_mni(ord,:);

    fprintf('\n  %s\n', label);
    fprintf('  %s\n', repmat('-', 1, numel(label)));
    fprintf('  Global max value   : %.6f\n', maxval);
    fprintf('  Tolerance (relative): %.1f%%  -> threshold = %.6f\n', ...
        tol_frac*100, thresh);
    fprintf('  Number of peaks    : %d\n', size(peak_pos,1));
    for k = 1:size(peak_pos,1)
        fprintf('    Peak %d:  native [%6.1f  %6.1f  %6.1f] mm  ->  MNI [%4.0f  %4.0f  %4.0f]  (val=%.6f)\n', ...
            k, peak_pos(k,1), peak_pos(k,2), peak_pos(k,3), ...
            peak_mni(k,1), peak_mni(k,2), peak_mni(k,3), peak_vals(k));
    end
end
