
%% beamforming in fieldtrip using BEM leadfields brain and spine

clear all
close all
clc

addpath('D:\brainspineconnectivity\source')
addpath('D:\spm')
spm('defaults','EEG')

sub='OP00212';
generic_dir = 'D:\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_realistic_cervical_brainspine_fine.mat');

HFC=1;
which_ori='all';


datwithEMGmerged = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);


load(geomfile)
D=spm_eeg_load(datwithEMGmerged);
grad_mm=D.sensors('MEG');
ftdat=spm2fieldtrip(D);

badchans=D.chanlabels(D.badchannels);

%remove bad channels here.
cfg=[];
cfg.channel=setdiff(ftdat.label,badchans);
ftdat=ft_selectdata(cfg,ftdat);

[Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'complex_bone_fields_cervical_brainspine_finegrid'));


sources_combined.pos    = [sources_cent.pos; sources_brain.pos];
sources_combined.inside = [logical(sources_cent.inside); sources_brain.inside];
sources_combined.unit   = sources_cent.unit;  % keep same unit


nsourcepoints = size(Gx,1);
nchannels     = size(Gx,2);

Lf.pos    = sources_combined.pos;                  % nsourcepoints x 3
Lf.inside = sources_combined.inside;     % all points inside
Lf.unit   = 'mm';
Lf.label  = grad_mm.label;              % nchannels x 1 cell
Lf.leadfielddimord = '{pos}_chan_ori';

% Allocate leadfield cell array
Lf.leadfield = cell(1,nsourcepoints);

for k = 1:nsourcepoints
    % Combine X/Y/Z components
    Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
end


cfg=[];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [10 35];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'no';
freqdat=ft_freqanalysis(cfg,ftdat);


cfg                     = [];
cfg.method              = 'infinite';
cfg.siunits=1;
cfg.grad=grad_mm;
cfg.conductivity = 1;

dummyvol = ft_prepare_headmodel(cfg,mesh_torso);


%construct filter for both conditions
cfg=[];
cfg.grid = sources_combined;
cfg.headmodel=dummyvol;
cfg.sourcemodel.leadfield=Lf;
cfg.dics.keepfilter='yes';
cfg.method = 'dics'; 
%cfg.fixedori='yes';
source_all = ft_sourceanalysis(cfg,freqdat); 


%% separate conditions

statidx=find(ftdat.trialinfo==1);
restidx=find(ftdat.trialinfo==2);

cfg=[];
cfg.trials=statidx;
statdat=ft_selectdata(cfg,ftdat);

cfg.trials=restidx;
restdat=ft_selectdata(cfg,ftdat);

cfg=[];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [10 35];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';
freqstat=ft_freqanalysis(cfg,statdat);
freqrest=ft_freqanalysis(cfg,restdat);

%now apply filter to each condition
cfg=[];
cfg.grid = sources_combined;
cfg.headmodel=dummyvol;
cfg.sourcemodel.leadfield=Lf;
cfg.dics.filter=source_all.avg.filter;
cfg.method = 'dics';
cfg.rawtrial='yes';

source_stat = ft_sourceanalysis(cfg,freqstat); %struct with single trials
source_rest = ft_sourceanalysis(cfg,freqrest);

cfg.rawtrial='no'; % struct with avg vals
source_stat_avg=ft_sourceanalysis(cfg,freqstat);
source_rest_avg=ft_sourceanalysis(cfg,freqrest);

n_trials=min([length(source_rest.trial) length(source_stat.trial)]);

sourcepow_trials_stat=nan(nsourcepoints, n_trials);
sourcepow_trials_rest=nan(nsourcepoints, n_trials);

for k=1:n_trials
    sourcepow_trials_stat(:,k)=source_stat.trial(k).pow;%sourcepoints x trials
    sourcepow_trials_rest(:,k)=source_rest.trial(k).pow;
end

%% permutation testing
 n_permutations = 100;

    tvals = zeros(1, nsourcepoints);
    pvals = zeros(1, nsourcepoints);

    for k = 1:nsourcepoints
        statdat = log(sourcepow_trials_stat(k,:));
        restdat = log(sourcepow_trials_rest(k,:));

        [~, p, ~, stats] = ttest(statdat, restdat);
        tvals(k) = stats.tstat;
        pvals(k) = p;
    end

    null_tvals = zeros(nsourcepoints, n_permutations);

    for i = 1:n_permutations
        for k = 1:nsourcepoints
            combined = [log(sourcepow_trials_stat(k,:)) log(sourcepow_trials_rest(k,:))];
            shuffled = combined(randperm(length(combined)));

            % Split back into two groups
            group1 = shuffled(1:n_trials);
            group2 = shuffled(n_trials+1:end);

            % Perform t-test
            [~, ~, ~, stats] = ttest(group1, group2);
            null_tvals(k, i) = stats.tstat;
        end
    end

    thresholds = prctile(null_tvals, 95, 2); % 1 for each sourcepoint
    thresholds_neg=prctile(null_tvals,5,2);
    %[pk, peakind] = max(tvals);

%%------------------------------------------------
% [~, ~, ~, stats] = ttest(sourcepow_trials_stat', sourcepow_trials_rest', 'Tail', 'right');
% 
% 
% source_diff.avg.pow=stats.tstat';

%--------------------------------------
%% to plot difference between mean across trials
source_diff=source_stat_avg; %copy
source_diff.avg.pow = log(source_stat_avg.avg.pow) - log(source_rest_avg.avg.pow);
%source_diff.avg.pow = (source_stat_avg.avg.pow - source_rest_avg.avg.pow) ./ source_rest_avg.avg.pow * 100;


%sep back to two source structures for plotting

n_spine = length(sources_cent.pos); 
n_total = length(source_diff.pos);   

spine_idx = 1:n_spine;
brain_idx = (n_spine+1):n_total;

%% --- Spine sources---
source_spine = struct();
source_spine.freq   = source_diff.freq;
source_spine.cfg    = source_diff.cfg;
source_spine.inside = source_diff.inside(spine_idx);
source_spine.pos    = source_diff.pos(spine_idx,:);
source_spine.unit   = source_diff.unit;
source_spine.method = source_diff.method;
%source_spine.avg.pow    = source_diff.avg.pow(spine_idx);  
source_spine.avg.pow=tvals(spine_idx);

%% --- Brain sources---
source_brain = struct();
source_brain.freq   = source_diff.freq;
source_brain.cfg    = source_diff.cfg;
source_brain.inside = source_diff.inside(brain_idx);
source_brain.pos    = source_diff.pos(brain_idx,:);
source_brain.unit   = source_diff.unit;
source_brain.method = source_diff.method;
%source_brain.avg.pow    = source_diff.avg.pow(brain_idx); 
source_brain.avg.pow=tvals(brain_idx);

% interpolate on meshes
cfg = [];
cfg.parameter = {'pow'};
brain_int = ft_sourceinterpolate(cfg, source_brain, mesh_brain); 
spine_int=ft_sourceinterpolate(cfg,source_spine, mesh_wm);

% brain_int.mask = (brain_int.pow > 0);% & (source_int.pow < source_int.thresh);
% spine_int.mask = (spine_int.pow > 0);% & (spine_int.pow < spine_int.thresh);

% [minval, minidx] = min(source_brain.avg.pow);
% minpos = source_diff.pos(minidx, :);
% 
% [~,idx] = min(vecnorm(mesh_brain.vertices - minpos,2,2));
% surfpos = mesh_brain.vertices(idx,:);
% 
% [minval, minidxsp] = min(source_spine.avg.pow);
% minpossp = source_spine.pos(minidxsp, :);
% 
% [~,idx] = min(vecnorm(mesh_wm.vertices - minpossp,2,2));
% surfpossp = mesh_wm.vertices(idx,:);

%% source plots


figure

%subplot(121)
cfg = [];
cfg.figure='gcf';
cfg.method = 'surface';
cfg.funparameter = 'pow';
% cfg.maskparameter = 'mask';
% cfg.maskstyle='opacity';
cfg.funcolormap = 'jet';%flipud(brewermap(512,'Blues')); 
%cfg.funcolorlim = 'zeromax';
cfg.projmethod = 'nearest';
cfg.surffile = mesh_brain;   

ft_sourceplot(cfg, brain_int); 
view(180,-8)
camlight



% hold on
% plot3(surfpos(1), surfpos(2), surfpos(3), 'ro', 'MarkerSize', 14, 'LineWidth', 3,...
%     'MarkerFaceColor','r')

%subplot(122)
cfg2=cfg;
%cfg2.funcolorlim = 'zeromax';
cfg2.projmethod = 'nearest';
cfg2.surffile = mesh_wm;   % your mesh struct with .pos and .tri
% cb2 = colorbar;
% cb2.Position = [0.85 0.15 0.02 0.3]; % manual placement
% ylabel(cb2,'t (spine)','FontSize',10)
ft_sourceplot(cfg2, spine_int);
view(80,25)
camlight
% hold on
% plot3(surfpossp(1), surfpossp(2), surfpossp(3), 'ro', 'MarkerSize', 12, 'LineWidth', 3,...
%     'MarkerFaceColor', 'r')



figure; plot(sources_cent.pos(:,2), tvals(spine_idx),'LineWidth',2)
hold on; plot(sources_cent.pos(:,2), thresholds(spine_idx),'k--')
hold on; plot(sources_cent.pos(:,2), thresholds_neg(spine_idx),'k--')
%%
figure; plot(sources_cent.pos(:,2), source_diff.avg.pow(spine_idx),'LineWidth',2)

%%
figure;

% Brain mesh: semitransparent red/pink
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.5, 'edgecolor', 'none');
hold on;

% Torso mesh: transparent blue shell
ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');

% White matter mesh: solid white/grey core
ft_plot_mesh(mesh_wm, 'facecolor', [0.9 0.9 0.9], 'facealpha', 1, 'edgecolor', 'none');

% Bone mesh: light beige, semitransparent
ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');

axis equal; 
camlight; lighting gouraud
%%
figure
ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');
hold on
plot3(sources_combined.pos(:,1), sources_combined.pos(:,2), sources_combined.pos(:,3),'b.')

