%% Define colour palette (RGB scaled 0–1)
colCoral  = [246 110  91] / 255;   % #F66E5B
colIndigo = [ 27  12  65] / 255;   % #1B0C41
colTeal   = [ 42 157 143] / 255;   % #2A9D8F
colOchre  = [201 162  39] / 255;   % #C9A227

%% 8) Plot
figure('Color','w','Position',[100 100 980 450]); hold on;

b = bar(binCenters, counts, 'grouped');

% Assign axis colours
b(1).FaceColor = colTeal;    % L–R
b(2).FaceColor = colCoral;   % C–C (dominant)
b(3).FaceColor = colOchre;   % D–V

% Optional: remove bar edges for cleaner vector look
for k = 1:3
    b(k).EdgeColor = 'none';
end

% 95th percentile threshold line
xline(thr95_max, '--', ...
    'Color', colIndigo, ...
    'LineWidth', 2);

% observed sig-max line
xline(obsMax, '-', ...
    'Color', colIndigo, ...
    'LineWidth', 2);

% ---- Highlight the observed bin (C–C bar only) ----
obsBin = discretize(obsMax, edges);

if ~isnan(obsBin)
    xObsBar = b(2).XEndPoints(obsBin);
    yObsBar = counts(obsBin,2);

    plot(xObsBar, yObsBar, 'o', ...
        'MarkerSize',10, ...
        'MarkerEdgeColor', colIndigo, ...
        'MarkerFaceColor', colCoral, ...
        'LineWidth',2);
end

xlabel('Permutation max coherence difference (binned)');
ylabel('Count (permutations)');

legend({'L–R dominant','C–C dominant','D–V dominant', ...
        '95th percentile (null max)', 'Observed sig-max'}, ...
       'Location','best', 'Box','off');

set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out');
box off;