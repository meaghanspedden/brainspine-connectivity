
%% source recon experiments in fieldtrip using BEM leadfields spine

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

%subs={'OP00212'};

subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

HFC=1;
LFop='spine';
correct_multcomp=1; 

subjResults=struct();

for ss=1:length(subs)
    sub=subs{ss};

    if strcmp(sub, 'OP00224')
    datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
    else
    datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D=spm_eeg_load(datwithEMGmerged);
    grad_mm=D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans=D.chanlabels(D.badchannels);

    %remove bad channels here.
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);

   spchanidx=find(grad_mm.coilpos(:,2) < 200);



    %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);

    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;     % all points inside
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);              % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';

    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end


    %% BEM and beamforming----------------------------
    %1. trial wise freq dat
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = [10 35];
    cfg.tapsmofrq  = 1;
    cfg.keeptrials='yes';
    freqdat_tr=ft_freqanalysis(cfg,ftdat); %trial wise freq


    %% 2. dummy head model for input config only
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);

    %% separate conditions

    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    cfg=[];
    cfg.trials=statidx;
    statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

    cfg.trials=restidx;
    restdat=ft_selectdata(cfg,freqdat_tr);

    %% get struct with mean freq data for constructing bf filter
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = [10 35];
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat=ft_freqanalysis(cfg,ftdat);

    %% construct filter using both conditions, use avg data
    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    %cfg.fixedori='yes';
    source_all = ft_sourceanalysis(cfg,freqdat);

    %% get trial-wise estimates with filter from above per condition
    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.filter=source_all.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.rawtrial='yes';
    cfg.keeptrials='yes';

    source_stat = ft_sourceanalysis(cfg,statdat); %struct per condition with single trials
    source_rest = ft_sourceanalysis(cfg,restdat);

    %% t test
    %extract to matrix
    [nTrials,~] = min([length(statidx) length(restidx)]);
    nSources     = size(source_stat.trial(1).pow,1);
    allpow_stat = zeros(nSources, nTrials);
    allpow_rest = zeros(nSources, nTrials);

    for t = 1:nTrials
        allpow_stat(:,t) = source_stat.trial(t).pow; % mean over freq
        allpow_rest(:,t) = source_rest.trial(t).pow;
    end

    mean_stat = mean(log(allpow_stat), 2); % mean over trials for plotting
    mean_rest = mean(log(allpow_rest), 2);

    %% plot mean over trials for each condition
%     figure
%     plot(sources_cent.pos(:,2),mean_rest,'LineWidth',1.5); hold on
%     plot(sources_cent.pos(:,2),mean_stat, 'LineWidth', 1.5)
%     legend({'rest', 'contraction'})
%     xlabel('Cranial caudal position')
%     title('BF and BEM mean each condition')

    %% ttest
    addpath 'C:\Users\mspedden\Documents\brainspineconnectivity'

    [H,P,CI,STATS] = ttest(log(allpow_stat'), log(allpow_rest'), 'tail', 'left');

    %correct p values
    [~, ~, ~, P]=fdr_bh(P,0.05,'pdep');

    sig_idx = P < 0.05;     % logical indices of significant points
    star_y = STATS.tstat(sig_idx) + 0.1*range(STATS.tstat);  % 10% above t-range

    x = sources_cent.pos(:,2);
    y = STATS.tstat;
    [~, I] = max(abs(y));

    figure('Color', 'w'); hold on; box on;
    plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 0.6]);
    plot(x(sig_idx), star_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(x(I), y(I), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    xlabel('Cranial-caudal distance (mm)', 'FontWeight', 'bold');
    ylabel('T statistic', 'FontWeight', 'bold');
    title('T statistics: rest vs contraction', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    xlim([min(x)-1, max(x)+1]);
    ylim([min(y)-0.1*range(y), max(y)+0.15*range(y)]);
    legend({'T values', 'p < 0.05', 'Max T'}, 'Location', 'best');
    set(gca, 'FontSize', 12);

    %% Interpolate on spine mesh
    source_spine = source_all;
    source_spine.avg.pow    = y'; %long

    % interpolate
    cfg = [];
    cfg.parameter = {'pow'};
    spine_int=ft_sourceinterpolate(cfg,source_spine, mesh_wm);


    %adjust torso mesh to plot less of lower body
    y=mesh_torso.vertices(:,2);
    keep_vert=y>-61;
    new_idx=zeros(size(keep_vert));
    new_idx(keep_vert)=1:sum(keep_vert);
    faces_keep=all(keep_vert(mesh_torso.faces),2);
    new_faces=new_idx(mesh_torso.faces(faces_keep,:));
    new_vertices=mesh_torso.vertices(keep_vert,:);
    mesh_cut.vertices=new_vertices;
    mesh_cut.faces=new_faces;
    mesh_cut.unit=mesh_torso.unit;

    % source plots
    % figure
    % cfg = [];
    % cfg.figure='gcf';
    % cfg.method = 'surface';
    % cfg.funparameter = 'pow';
    % cfg.funcolormap = 'magma';
    % cfg.funcolorlim='minzero';
    % %cfg.funcolorlim = [-.034 -.03];
    % cfg.projmethod = 'nearest';
    % cfg.surffile = mesh_wm;
    % ft_sourceplot(cfg, spine_int);
    % hold on
    % ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.2, 'edgecolor', 'none');
    % ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.07, 'edgecolor', 'none');
    % ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.1, 'edgecolor', 'none');
    % view( -250, -1)
    % camlight
    %% add a mask
    source_mask=source_spine; %copy
    source_mask.avg.pow = (P < 0.05)';

    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'nearest';
    source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_wm);
    spine_int.mask=source_mask_int.pow;

%     figure
%     cfg = [];
%     cfg.figure='gcf';
%     cfg.method = 'surface';
%     cfg.funparameter = 'pow';
%     cfg.funcolormap = 'parula';
%     cfg.maskparameter = 'mask';
%     cfg.maskstyle='opacity';
%     cfg.projmethod = 'nearest';
%     cfg.surffile = mesh_wm;
%     ft_sourceplot(cfg, spine_int);
%     camlight
%     hold on
%     ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.2, 'edgecolor', 'none');
%     ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.07, 'edgecolor', 'none');
%     ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.1, 'edgecolor', 'none');
%     view( -250, -1)
%     camlight
    %% permutation testing
    for t = 1:nTrials
        allpow_stat(:,t) = source_stat.trial(t).pow; % mean over freq
        allpow_rest(:,t) = source_rest.trial(t).pow;
    end
    %% permutation testing
    n_permutations = 100;

    null_tvals = zeros(nSources, n_permutations);

    for i = 1:n_permutations
        for k = 1:nSources
            combined = [log(allpow_stat(k,:)) log(allpow_rest(k,:))];
            shuffled = combined(randperm(length(combined)));

            % Split back into two groups
            group1 = shuffled(1:nTrials);
            group2 = shuffled(nTrials+1:end);

            % Perform t-test
            [~, ~, ~, stats] = ttest(group1, group2, 'tail', 'left');
            null_tvals(k, i) = stats.tstat;
            null_min_t(i) = min(stats.tstat); % store most extreme (minimum) t

        end
    end
    if correct_multcomp
        threshold_neg=prctile(null_min_t,5,2);
    else
        threshold_neg = prctile(null_tvals, 5,2);  % 5th percentile for negative tail
    end
%% Plot permutation results
f = figure('Color', 'w'); 
hold on; box on; grid on;

% === Permutation thresholds (shaded null region) ===
fill([x; flipud(x)], ...
     [threshold_neg.*ones(size(x)); flipud(zeros(size(x)))], ... % shaded area from threshold to 0
     [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

% Threshold line
plot(x, threshold_neg.*ones(size(x)), 'k--', 'LineWidth', 1);

% T-values
y = STATS.tstat';
plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 0.6]); 

% Significant points
plot(x(sig_idx), y(sig_idx), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);

xlabel('Cranial–caudal distance (mm)', 'FontWeight', 'bold');
ylabel('T statistic', 'FontWeight', 'bold');
title('Rest vs Contraction (left-tailed)', 'FontSize', 14, 'FontWeight', 'bold');

legend({'Permutation 5% null region', 'Threshold', 'Parametric T values', 'FDR corr p < 0.05'}, ...
       'Location', 'best');

xlim([min(x)-1, max(x)+1]);
ylim([min([y; threshold_neg])-0.1*range(y), max([y; 0])+0.15*range(y)]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

waitfor(f);

% Save results
subjResults(ss).tstat = STATS.tstat;
subjResults(ss).sig_idx = sig_idx;
subjResults(ss).thr95_neg = threshold_neg;
subjResults(ss).allpow_stat = allpow_stat;
subjResults(ss).allpow_rest = allpow_rest;

end
%% prevalence of sig diffs


nSubjects = length(subjResults);

% Initialize prevalence counters for left-tailed test
sig_param_neg = false(nSubjects,1);  % parametric t < 0
sig_perm_neg  = false(nSubjects,1);  % permutation t < threshold

for ss = 1:nSubjects
    t = subjResults(ss).tstat';
    sig_idx = subjResults(ss).sig_idx;
    thr_neg = subjResults(ss).thr95_neg;
    
    % Parametric significance (left-tailed)
    if any(sig_idx)
        sig_param_neg(ss) = true;
    end
    
    % Permutation significance (left-tailed)
    if any(t < thr_neg)
        sig_perm_neg(ss) = true;
    end
end

fprintf('Parametric (decreased effect): %d/%d subjects\n', sum(sig_param_neg), nSubjects);
fprintf('Permutation (decreased effect): %d/%d subjects\n', sum(sig_perm_neg), nSubjects);

%% visualise over subjects

% significance vectors
nSubjects = length(subjResults);
nSources = length(subjResults(1).tstat);
sigMat_param = false(nSubjects, nSources);
sigMat_perm  = false(nSubjects, nSources);

for ss = 1:nSubjects
    sigMat_param(ss,:) = subjResults(ss).sig_idx;
    sigMat_perm(ss,:)  = subjResults(ss).tstat' < subjResults(ss).thr95_neg;
end

% Compute prevalence
prevalence_param = mean(sigMat_param,1) * 100;  
prevalence_perm  = mean(sigMat_perm,1)  * 100;

nSubj = numel(subjResults);
x = sources_cent.pos(:,2);

figure; hold on;
cmap = [
    27,158,119;
    217,95,2;
    117,112,179;
    231,41,138;
    102,166,30;
    230,171,2;
    166,118,29;
    102,102,102;
    55,126,184
] / 255;

for s = 1:nSubj
    tvals = subjResults(s).tstat;
    thr   = subjResults(s).thr95_neg;  % permutation threshold (left-tailed)
    sig   = tvals' < thr;               % significant decreases
    
    % Choose color based on significance
    if any(sig)
        c = cmap(s,:);
    else
        c = [0.7 0.7 0.7];
    end
    
    % Plot segments — color only where both adjacent points are significant
    for i = 1:length(x)-1
        if sig(i) && sig(i+1)
            plot(x(i:i+1), tvals(i:i+1), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        else
            plot(x(i:i+1), tvals(i:i+1), '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
    end

    plot(x(sig), tvals(sig), '.', 'Color', c, 'MarkerSize', 12, 'HandleVisibility', 'off'); %sig pts
    h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 1.5);
end

yline(0, ':k', 'HandleVisibility', 'off');
xlabel('Cranial-caudal position (mm)');
ylabel('T-statistic');
title('Significant decreases (Static contraction < rest)');
legend(h, arrayfun(@(s) sprintf('Subj %d', s), 1:nSubj, 'UniformOutput', false), ...
       'Location', 'bestoutside');

set(gca, 'FontSize', 13)
grid on;



%% sorted by height 

heighttable=readtable('C:\Users\mspedden\Documents\SC_subs_heights.csv');
heights=heighttable.Var2;

[sortedHeights, sortIdx] = sort(heights, 'descend');  % tallest first
subjResultsSorted = subjResults(sortIdx);

cmapSorted = cmap(sortIdx, :);

nSubj = numel(subjResultsSorted);
figure; hold on;

for s = 1:nSubj
    tvals = subjResultsSorted(s).tstat;
    thr   = subjResultsSorted(s).thr95_neg;
    sig   = tvals' < thr;

    % choose color
    if any(sig)
        c = cmapSorted(s,:);
    else
        c = [0.7 0.7 0.7];
    end

    % plot segments
    for i = 1:length(x)-1
        if sig(i) && sig(i+1)
            plot(x(i:i+1), tvals(i:i+1), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        else
            plot(x(i:i+1), tvals(i:i+1), '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
    end

    plot(x(sig), tvals(sig), '.', 'Color', c, 'MarkerSize', 12, 'HandleVisibility', 'off'); % sig pts
    h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 1.5);  % for legend
end

yline(0, ':k', 'HandleVisibility', 'off');
xlabel('Cranial-caudal position (mm)');
ylabel('T-statistic');
title('Significant decreases sorted by height');
legend(h, arrayfun(@(s) sprintf('Subj %d', s), 1:nSubj, 'UniformOutput', false), 'Location', 'bestoutside');
set(gca, 'FontSize', 13)
grid on;


