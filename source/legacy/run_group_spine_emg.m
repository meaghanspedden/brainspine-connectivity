%% run_group_spine_only.m
clear all; close all; clc;

fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';

out_suffix = '_true_bem';
fwhm_mm    = 20;

addpath(spm_path);
addpath(fieldtrip_path);
ft_defaults;

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end



%% Load geometry
geom_exp   = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_bone    = geom_exp.mesh_bone;
mesh_lungs   = geom_exp.mesh_lungs;
mesh_heart   = geom_exp.mesh_heart;

geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;

roi_sensory = 14:20;
roi_motor   = 21:27;
roi_all     = 14:27;
z_offset    = 10;
cord_pos    = sources_cent.pos(:,2);

%% Colourmap
ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% Load group results
load(fullfile(save_dir, ['groupRes_spine_DICS_bemv2' out_suffix '.mat']));


run_group_spine(subjResults, sources_cent, cord_pos, save_dir, out_suffix, cmap_hot, ...
    mesh_wm, mesh_brain, mesh_torso, mesh_bone, mesh_lungs, mesh_heart, ...
    fwhm_mm, fig_dir, roi_sensory, roi_motor, roi_all, z_offset);

fprintf('Done.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function run_group_spine(subjResults, sources_cent, cord_pos, save_dir, out_suffix, cmap_hot, ...
                          mesh_wm, mesh_brain, mesh_torso, mesh_bone, mesh_lungs, mesh_heart, ...
                          fwhm, fig_dir, roi_sensory, roi_motor, roi_all, z_offset)

nSubjects     = length(subjResults);
nsourcepoints = size(sources_cent.pos,1);
all_masks     = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    m = double(subjResults(s).coh_diff > subjResults(s).thr95);
    all_masks(:,s) = m(:);
end

sig_pos = false(nSubjects,1);
for s = 1:nSubjects
    sig_pos(s) = any(all_masks(:,s));
end
fprintf('  %d/%d subjects show significant spine-EMG coherence (smoothed %d mm, true BEM)\n', ...
    sum(sig_pos), nSubjects, fwhm);

threshold      = 0.2;
prevalence_loc = mean(all_masks, 2);

group_ft         = [];
group_ft.pos     = sources_cent.pos;
group_ft.inside  = sources_cent.inside;
group_ft.pow     = prevalence_loc;
group_ft.pow(group_ft.pow < threshold) = 0;

cfg = []; cfg.parameter = 'pow'; cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);

mesh_cut = clip_torso(mesh_torso);

% Group prevalence — full figure with torso
hfig_prev = figure;
cfg2 = []; cfg2.method = 'surface'; cfg2.funparameter = 'pow';
cfg2.maskparameter = 'mask';
cfg2.funcolorlim   = [threshold max(group_int.pow)];
cfg2.funcolormap   = cmap_hot;
cfg2.projmethod    = 'nearest';
cfg2.surffile      = mesh_wm;
cfg2.opacitylim    = [threshold max(group_int.pow)];
cfg2.opacitymap    = 'rampup';
ft_sourceplot(cfg2, group_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14; hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.3, 'edgecolor','none');
ft_plot_mesh(mesh_lungs,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_heart,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
% ROI markers — sensory (blue) and motor (red), offset in z to sit on surface
for s = roi_sensory
    plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
        sources_cent.pos(s,3) + z_offset, ...
        'o', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 10, 'LineWidth', 1);
end
for s = roi_motor
    plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
        sources_cent.pos(s,3) + z_offset, ...
        'o', 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 10, 'LineWidth', 1);
end
% dummy handles for legend
plot3(nan,nan,nan,'o','MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','w',...
    'MarkerSize',10,'DisplayName','Sensory (C6-C8)');
plot3(nan,nan,nan,'o','MarkerFaceColor',[0.8 0.2 0.2],'MarkerEdgeColor','w',...
    'MarkerSize',10,'DisplayName','Motor (C8-T1)');
legend('Location','best','FontSize',10);
title(sprintf('Group prevalence — spine-EMG (smoothed %d mm, true BEM)', fwhm), ...
    'Interpreter','none');
savefig(hfig_prev, fullfile(fig_dir, ['step2_group_spineEMG_prevalence' out_suffix '.fig']));
%close(hfig_prev);

% Group prevalence — mesh only with star at peak
hfig_grp_mesh = figure('Color','w','Position',[100 100 400 650]);
cfg3 = []; cfg3.method = 'surface'; cfg3.funparameter = 'pow';
cfg3.funcolorlim = [threshold max(group_int.pow)];
cfg3.funcolormap = cmap_hot;
cfg3.projmethod  = 'nearest';
cfg3.surffile    = mesh_wm;
cfg3.opacitylim  = [threshold max(group_int.pow)];
cfg3.opacitymap  = 'rampup';
ft_sourceplot(cfg3, group_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
title(sprintf('Group prevalence — spine-EMG (smoothed %d mm, true BEM)', fwhm), ...
    'Interpreter','none','FontSize',13);
prev_vals = group_ft.pow; prev_vals(prev_vals < threshold) = -inf;
if any(isfinite(prev_vals))
    peaks_grp = sources_cent.pos(prev_vals >= ...
        max(prev_vals(isfinite(prev_vals)))*0.99, :);
    hold on;
    scatter3(peaks_grp(:,1), peaks_grp(:,2), peaks_grp(:,3)+10, ...
        200, 'p', 'filled', ...
        'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k','LineWidth',1.5);
    scatter_obj = findobj(gca,'Type','Scatter');
    uistack(scatter_obj,'top');
end
hold on
% ROI markers — sensory (blue) and motor (red), offset in z to sit on surface
for s = roi_sensory
    plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
        sources_cent.pos(s,3) + z_offset, ...
        'o', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 10, 'LineWidth', 1);
end
for s = roi_motor
    plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
        sources_cent.pos(s,3) + z_offset, ...
        'o', 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 10, 'LineWidth', 1);
end
% dummy handles for legend
plot3(nan,nan,nan,'o','MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','w',...
    'MarkerSize',10,'DisplayName','Sensory (C6-C8)');
plot3(nan,nan,nan,'o','MarkerFaceColor',[0.8 0.2 0.2],'MarkerEdgeColor','w',...
    'MarkerSize',10,'DisplayName','Motor (C8-T1)');
legend('Location','best','FontSize',10);
savefig(hfig_grp_mesh, fullfile(fig_dir, ...
    ['step2_group_spineEMG_prevalence_meshonly' out_suffix '.fig']));
%close(hfig_grp_mesh);

% Subject line plot
subj_cmap = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; ...
             230,171,2;  166,118,29; 102,102,102; 55,126,184] / 255;
x = sources_cent.pos(:,2);
hfig_lines = figure; hold on;
for s = 1:nSubjects
    cdiff = subjResults(s).coh_diff;
    thr   = subjResults(s).thr95;
    sig   = cdiff > thr;
    c     = subj_cmap(s,:);
    for i = 1:length(x)-1
        if sig(i) && sig(i+1)
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', c, ...
                'LineWidth',1.5,'HandleVisibility','off');
        else
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color',[0.7 0.7 0.7], ...
                'LineWidth',1,'HandleVisibility','off');
        end
    end
    plot(x(sig), cdiff(sig), '.', 'Color', c, 'MarkerSize',12,'HandleVisibility','off');
    if sig_pos(s)
        h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth',1.5);
    else
        h(s) = plot(nan, nan, '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.5);
    end
end
% ROI shading on line plot
% ROI shading — y limits from all subjects combined
all_coh = cat(1, subjResults(:).coh_diff);
yl_pad  = [min(all_coh)*1.1, max(all_coh)*1.1];
fill([cord_pos(roi_all(1))   cord_pos(roi_all(end)) ...
      cord_pos(roi_all(end)) cord_pos(roi_all(1))], ...
     [yl_pad(1) yl_pad(1) yl_pad(2) yl_pad(2)], ...
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha', 0.3,'DisplayName','ROI (C6-T1)');

yline(0,':k','HandleVisibility','off');
xlabel('Cranial-caudal position (mm)'); ylabel('Coherence difference (stat-rest)');
title(sprintf('Spine-EMG coherence differences (smoothed %d mm, true BEM)', fwhm));
legend(h, arrayfun(@(s) sprintf('Participant %d',s), 1:nSubjects,'UniformOutput',false), ...
    'Location','bestoutside');
set(gca,'FontSize',13); grid on;
savefig(hfig_lines, fullfile(fig_dir, ...
    ['step2_group_spineEMG_subject_lines' out_suffix '.fig']));
%close(hfig_lines);

% Cluster ROI for VE
mask_thresh = prevalence_loc >= threshold;
pos_thresh  = sources_cent.pos(mask_thresh,:);
if ~isempty(pos_thresh)
    distMat     = squareform(pdist(pos_thresh));
    G           = graph(distMat < 6);
    bins        = conncomp(G);
    [~, idxMax] = max(histcounts(bins, 1:(max(bins)+1)));
    ROIpos      = pos_thresh(bins == idxMax, :);

    hfig_clust = figure; hold on;

% ROI shading
yl_clust = [0, max(prevalence_loc)*1.1];
fill([cord_pos(roi_all(1))   cord_pos(roi_all(end)) ...
      cord_pos(roi_all(end)) cord_pos(roi_all(1))], ...
     [yl_clust(1) yl_clust(1) yl_clust(2) yl_clust(2)], ...
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3, ...
     'HandleVisibility','off');

plot(sources_cent.pos(:,2), prevalence_loc);
for k = 1:size(ROIpos,1)
    plot(ROIpos(k,2), 0.2, 'r*');
end

    hfig_clust = figure; hold on;
    plot(sources_cent.pos(:,2), prevalence_loc);
    for k = 1:size(ROIpos,1)
        plot(ROIpos(k,2), 0.2, 'r*');
    end
    xlabel('Cranial-caudal position (mm)'); ylabel('Prevalence');
    title('Spinal prevalence + VE cluster (true BEM)','Interpreter','none');
    savefig(hfig_clust, fullfile(fig_dir, ...
        ['step2_spinal_prevalence_VE_cluster' out_suffix '.fig']));
    %close(hfig_clust);

    save(fullfile(save_dir, ['cluster_spineEMG_pos_bemv2' out_suffix '.mat']), 'ROIpos');
    fprintf('  ROI cluster saved: %d sources\n', size(ROIpos,1));
else
    fprintf('  WARNING: no sources exceeded prevalence threshold — ROI cluster not saved\n');
end

end

% -------------------------------------------------------------------------
function [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm)
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
end

% -------------------------------------------------------------------------
function thr = compute_threshold(cohDiff_perm, mult_comp_corr, nsourcepoints)
    maxPerm = max(cohDiff_perm, [], 1);
    if mult_comp_corr
        thr = prctile(maxPerm, 95);
    else
        thr = prctile(cohDiff_perm, 95, 2);
        warning('Uncorrected per-source threshold used.');
    end
end

% -------------------------------------------------------------------------
function pvals = compute_pvals(coh_diff, cohDiff_perm)
    nPerm = size(cohDiff_perm, 2);
    pvals = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (nPerm + 1);
end

% -------------------------------------------------------------------------
function invp_smooth = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm)
    invp_smooth = zeros(nsourcepoints,1);
    for s = 1:nsourcepoints
        permDist = sort(cohDiff_perm(s,:));
        obsVal   = coh_diff(s);
        xgrid    = linspace(min(permDist), max(permDist), 200);
        p_emp    = arrayfun(@(x) (sum(permDist >= x)+1)/(nPerm+1), xgrid);
        logp_sm  = smooth(xgrid, -log10(p_emp), 0.15, 'loess');
        obsVal_c = min(max(obsVal, xgrid(1)), xgrid(end));
        invp_smooth(s) = interp1(xgrid, logp_sm, obsVal_c, 'linear');
    end
end

% -------------------------------------------------------------------------
function mesh_cut = clip_torso(mesh_torso)
    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx   = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep        = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;
end

% -------------------------------------------------------------------------
function W = make_gaussian_smoother(pos_mm, fwhm_mm, radius_mm)
    sigma = fwhm_mm / 2.355;
    if nargin < 3 || isempty(radius_mm), radius_mm = 3*sigma; end
    N   = size(pos_mm, 1);
    Mdl = KDTreeSearcher(pos_mm);
    [idx, dist] = rangesearch(Mdl, pos_mm, radius_mm);
    ii = []; jj = []; vv = [];
    for i = 1:N
        j = idx{i}; d = dist{i};
        w = exp(-0.5*(d./sigma).^2);
        ii = [ii; repmat(i,numel(j),1)]; jj = [jj; j(:)]; vv = [vv; w(:)];
    end
    W  = sparse(ii,jj,vv,N,N);
    rs = full(sum(W,2)); rs(rs==0) = 1;
    W  = spdiags(1./rs,0,N,N) * W;
end