%% QUANTITATIVE COMPARISON OF THREE LEADFIELDS — corrected channel ordering
% LF1: bem_v1  (leadfield_cord)
% LF2: bem_v2  (leadfield_cord)
% LF3: build_leadfield_matrices (Gx/Gy/Gz -> Lf)

clear all; close all; clc;

addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')
addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile    = fullfile(generic_dir, 'geometries_cervical_realistic.mat');
LFop        = 'spine';

load(geomfile)

% --- shared channel index (grad_mm convention = LF3 convention) ---
sub              = 'OP00212';
datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
    ['sub-' sub], 'ses-001', 'meg', ...
    'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
D       = spm_eeg_load(datwithEMGmerged);
grad_mm = D.sensors('MEG');
spchanidx = find(grad_mm.coilpos(:,2) < 200);   % spinal channels in grad_mm

%% ---- Load LF1: BEM v1 + build position-based channel map ----
load('C:\Users\mspedden\Documents\bem_spine_fields\bem_v1_leadfield_cervical_realistic_bem_.mat')
grad_v1 = leadfield_cord.cfg.grad;

% For each spinal channel in grad_mm, find matching row in grad_v1 (metres->mm)
spc_in_v1 = zeros(numel(spchanidx), 1);
for i = 1:numel(spchanidx)
    p  = grad_mm.coilpos(spchanidx(i), :);           % mm
    d  = sqrt(sum((grad_v1.coilpos*1000 - p).^2, 2));
    [~, spc_in_v1(i)] = min(d);
end

LF1 = leadfield_cord;
LF1.label = grad_mm.label(spchanidx);   % adopt grad_mm convention
for i = 1:numel(leadfield_cord.leadfield)
    if ~isempty(leadfield_cord.leadfield{i})
        LF1.leadfield{i} = leadfield_cord.leadfield{i}(spc_in_v1, :);
    end
end

%% ---- Load LF2: BEM v2 + build position-based channel map ----
load('C:\Users\mspedden\Documents\bem_spine_fields\bem_v2_leadfield_cervical_realistic_bem_.mat')
grad_v2 = leadfield_cord.cfg.grad;

spc_in_v2 = zeros(numel(spchanidx), 1);
for i = 1:numel(spchanidx)
    p  = grad_mm.coilpos(spchanidx(i), :);
    d  = sqrt(sum((grad_v2.coilpos*1000 - p).^2, 2));
    [~, spc_in_v2(i)] = min(d);
end

LF2 = leadfield_cord;
LF2.label = grad_mm.label(spchanidx);
for i = 1:numel(leadfield_cord.leadfield)
    if ~isempty(leadfield_cord.leadfield{i})
        LF2.leadfield{i} = leadfield_cord.leadfield{i}(spc_in_v2, :);
    end
end

%% ---- Load LF3: Gx/Gy/Gz (already in grad_mm convention) ----
[Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);
nsourcepoints = size(Gx, 1);
Gx = Gx(:, spchanidx);
Gy = Gy(:, spchanidx);
Gz = Gz(:, spchanidx);

LF3.pos             = sources_cent.pos;
LF3.inside          = sources_cent.inside;
LF3.unit            = 'mm';
LF3.label           = grad_mm.label(spchanidx);
LF3.leadfielddimord = '{pos}_chan_ori';
LF3.leadfield       = cell(1, nsourcepoints);
for k = 1:nsourcepoints
    LF3.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)'];
end

%% ---- Sanity check: channel labels consistent across all three ----
assert(isequal(LF1.label, LF2.label), 'LF1 and LF2 labels differ!');
assert(isequal(LF1.label, LF3.label), 'LF1 and LF3 labels differ!');
fprintf('Channel order check passed — all three LFs use same label order.\n');

inside_idx = find(sources_cent.inside);
nInside    = numel(inside_idx);
%% ---- Sanity check: print first 5 channels + one source values ----
s_test = find(sources_cent.inside, 1, 'first');
s_test = inside_idx(round(numel(find(sources_cent.inside))/2));

fprintf('\n%-10s %14s %14s %14s\n', 'channel', 'LF1 (v1)', 'LF2 (v2)', 'LF3 (Gxyz)');
for ch = 1:5
    fprintf('%-10s %14.4e %14.4e %14.4e\n', ...
        grad_mm.label{spchanidx(ch)}, ...
        LF1.leadfield{s_test}(ch,1), ...
        LF2.leadfield{s_test}(ch,1), ...
        LF3.leadfield{s_test}(ch,1));
end

%% ---- Inside indices ----

lf_names   = {'BEM v1', 'BEM v2', 'Gx/Gy/Gz'};
LFs        = {LF1, LF2, LF3};
ypos       = sources_cent.pos(inside_idx, 2);

%% ---- Per-source Frobenius norm and condition number ----
frob     = zeros(nInside, 3);
cond_num = zeros(nInside, 3);

for lv = 1:3
    for si = 1:nInside
        s  = inside_idx(si);
        lf = LFs{lv}.leadfield{s};
        frob(si,lv)     = norm(lf, 'fro');
        cond_num(si,lv) = cond(lf);
    end
end

%% ---- Pairwise correlations using normalised leadfields ----
% Normalise each LF by its mean Frobenius norm before correlating
% so that unit differences don't drive the correlation
LFs_norm = LFs;
for lv = 1:3
    scale = mean(frob(:,lv));
    for si = 1:nInside
        s = inside_idx(si);
        LFs_norm{lv}.leadfield{s} = LFs{lv}.leadfield{s} / scale;
    end
end

pairs      = [1 2; 1 3; 2 3];
pair_names = {'v1 vs v2', 'v1 vs Gxyz', 'v2 vs Gxyz'};
pair_cols  = [0.5 0.1 0.8; 0.9 0.6 0.1; 0.1 0.6 0.6];

r_mat = zeros(nInside, 3);

fprintf('\n========== LEADFIELD SUMMARY ==========\n');
fprintf('%-25s %12s %12s %12s\n', 'Metric', lf_names{:});
fprintf('%s\n', repmat('-',1,65));
metrics = {'Mean Frob norm','Std Frob norm','Mean cond number','Median cond number'};
vals    = [mean(frob); std(frob); mean(cond_num); median(cond_num)];
for m = 1:size(vals,1)
    fprintf('%-25s %12.4e %12.4e %12.4e\n', metrics{m}, vals(m,1), vals(m,2), vals(m,3));
end

fprintf('\n========== PAIRWISE CORRELATIONS (normalised) ==========\n');
for pp = 1:3
    a = pairs(pp,1); b = pairs(pp,2);
    for si = 1:nInside
        s   = inside_idx(si);
        lfa = LFs_norm{a}.leadfield{s}(:);
        lfb = LFs_norm{b}.leadfield{s}(:);
        r_mat(si,pp) = corr(lfa, lfb);
    end
    fprintf('%-15s  mean r = %.4f,  median r = %.4f,  min r = %.4f\n', ...
        pair_names{pp}, mean(r_mat(:,pp)), median(r_mat(:,pp)), min(r_mat(:,pp)));
end

%% ---- Figure 1: Frobenius norm along cord ----
cols = [0.2 0.5 0.9; 0.9 0.3 0.2; 0.2 0.7 0.3];

figure('Color','w','Position',[100 100 900 350]); hold on;
for lv = 1:3
    plot(ypos, frob(:,lv), '-', 'Color', cols(lv,:), 'LineWidth', 1.5, ...
        'DisplayName', lf_names{lv});
end
xlabel('Cranio–caudal position (mm)', 'FontSize',13);
ylabel('Frobenius norm', 'FontSize',13);
title('Leadfield magnitude along spinal cord', 'FontSize',13);
legend('Location','best','FontSize',12); grid on; box off; set(gca,'FontSize',13);

%% ---- Figure 2: Condition number along cord ----
figure('Color','w','Position',[100 100 900 350]); hold on;
for lv = 1:3
    plot(ypos, log10(cond_num(:,lv)), '-', 'Color', cols(lv,:), 'LineWidth', 1.5, ...
        'DisplayName', lf_names{lv});
end
xlabel('Cranio–caudal position (mm)', 'FontSize',13);
ylabel('log_{10}(condition number)', 'FontSize',13);
title('Leadfield condition number along spinal cord', 'FontSize',13);
legend('Location','best','FontSize',12); grid on; box off; set(gca,'FontSize',13);

%% ---- Figure 3: Pairwise normalised difference along cord ----
figure('Color','w','Position',[100 100 900 350]); hold on;
for pp = 1:3
    a = pairs(pp,1); b = pairs(pp,2);
    rel_diff = zeros(nInside,1);
    for si = 1:nInside
        s   = inside_idx(si);
        lfa = LFs_norm{a}.leadfield{s};
        lfb = LFs_norm{b}.leadfield{s};
        rel_diff(si) = norm(lfa-lfb,'fro') / (0.5*(norm(lfa,'fro')+norm(lfb,'fro')));
    end
    plot(ypos, rel_diff, '-', 'Color', pair_cols(pp,:), 'LineWidth', 1.5, ...
        'DisplayName', pair_names{pp});
end
xlabel('Cranio–caudal position (mm)', 'FontSize',13);
ylabel('Normalised difference', 'FontSize',13);
title('Pairwise leadfield similarity along spinal cord (normalised)', 'FontSize',13);
legend('Location','best','FontSize',12); grid on; box off; set(gca,'FontSize',13);

%% ---- Figure 4: Pairwise correlation along cord ----
figure('Color','w','Position',[100 500 900 350]); hold on;
for pp = 1:3
    plot(ypos, r_mat(:,pp), '-', 'Color', pair_cols(pp,:), 'LineWidth', 1.5, ...
        'DisplayName', pair_names{pp});
end
yline(0, ':k', 'HandleVisibility','off');
xlabel('Cranio–caudal position (mm)', 'FontSize',13);
ylabel('Correlation (r)', 'FontSize',13);
title('Pairwise leadfield correlation along spinal cord (normalised)', 'FontSize',13);
legend('Location','best','FontSize',12); grid on; box off; set(gca,'FontSize',13);

%% ---- Figure 5: Boxplot of raw Frobenius norms ----
figure('Color','w','Position',[100 500 500 400]);
boxplot(frob, lf_names, 'Symbol','', 'Widths',0.5);
ylabel('Frobenius norm', 'FontSize',13);
title('Distribution of LF magnitudes (inside sources)', 'FontSize',13);
set(gca,'FontSize',13); grid on; box off;