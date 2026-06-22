%% =========================================================================
%  FIGURE: Overlaid aligned spectra
%% =========================================================================
rel_win  = 10;
rel_res  = 0.25;
rel_axis = -rel_win : rel_res : rel_win;
smooth_w = 5;

pair_cols    = {col_BE, col_BS, col_SE};
pair_names   = {'Brain-EMG', 'Brain-Spine', 'Spine-EMG'};
band_data    = {'be_band', 'bs_band', 'se_band'};
match_fields = {[], 'match_BS', 'match_SE'};
match_hz     = {[], 'bs_match_hz', 'se_match_hz'};

yl_patch = [-2 4];

hfig_ov = figure('Color','w','Position',[100 100 1100 380]);

for pp = 1:3
    ax = subplot(1,3,pp); hold on;
    col = pair_cols{pp};

    mat = nan(nSubs, numel(rel_axis));

    for ss = 1:nSubs
        r        = results(ss);
        spec_raw = r.(band_data{pp});
        spec_sm  = movmean(spec_raw, smooth_w);
        spec_norm = zscore(spec_sm);

        fax_rel  = r.freqs_band - r.be_peak_hz;
        in_range = rel_axis >= min(fax_rel) & rel_axis <= max(fax_rel);
        if sum(in_range) < 2, continue; end
        mat(ss, in_range) = interp1(fax_rel, spec_norm, rel_axis(in_range), 'linear');
    end

    % Shaded ±2 Hz window
    patch([-match_window_hz match_window_hz match_window_hz -match_window_hz], ...
        [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
        [0.85 0.85 0.85], 'FaceAlpha', 0.4, 'EdgeColor','none', ...
        'HandleVisibility','off');

    % Individual traces
    for ss = 1:nSubs
        y = mat(ss,:);
        if all(isnan(y)), continue; end
        plot(rel_axis, y, '-', 'Color', [col 0.20], ...
            'LineWidth', 0.8, 'HandleVisibility','off');
    end

    % Group median
    med_trace = median(mat, 1, 'omitnan');
    ok        = ~isnan(med_trace);
    plot(rel_axis(ok), med_trace(ok), '-', 'Color', col, ...
        'LineWidth', 3.0, 'HandleVisibility','off');

    % Mark matched peaks
    for ss = 1:nSubs
        r = results(ss);
        if pp == 1
            [~, i0] = min(abs(rel_axis));
            plot(0, mat(ss, i0), '.', 'Color', col, ...
                'MarkerSize', 14, 'HandleVisibility','off');
        else
            mf = match_fields{pp};
            mh = match_hz{pp};
            if r.(mf)
                rel_pos   = r.(mh) - r.be_peak_hz;
                [~, irel] = min(abs(rel_axis - rel_pos));
                y_val     = mat(ss, irel);
                if ~isnan(y_val)
                    scatter(rel_pos, y_val, 60, col, 'o', 'filled', ...
                        'MarkerEdgeColor','k', 'LineWidth',0.8, ...
                        'HandleVisibility','off');
                end
            else
                [~, i0] = min(abs(rel_axis));
                y_val   = mat(ss, i0);
                if ~isnan(y_val)
                    scatter(0, y_val, 40, [0.6 0.6 0.6], 'x', ...
                        'LineWidth', 1.5, 'HandleVisibility','off');
                end
            end
        end
    end

    xline(0, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    xlim([-rel_win rel_win]);
    ylim(yl_patch);
    xlabel('Hz relative to brain-EMG peak', 'FontSize', 13);
    if pp == 1
        ylabel('Coherence (z-score)', 'FontSize', 13);
    end
    title(pair_names{pp}, 'FontSize', 13, 'FontWeight','normal');

    if pp > 1
        mf      = match_fields{pp};
        n_match = sum([results.(mf)]);
        text(rel_win*0.95, 3.7, sprintf('%d/%d matched', n_match, nSubs), ...
            'HorizontalAlignment','right', 'FontSize', 11, 'Color', col);
    end

    set(gca, 'FontSize', 13); box off;
end

sgtitle('Coherence spectra aligned to brain-EMG peak frequency', 'FontSize', 13);
print(hfig_ov, fullfile(fig_dir, 'overlaid_aligned_spectra.png'), '-dpng', '-r300');
fprintf('Figure saved.\n');