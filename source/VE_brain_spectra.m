% virtual electrode brain(from prevalence analysis cluster) 
clear all
close all
clc

addpath('C:\Users\mspaedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end


subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
 'OP00225', 'OP00221', 'OP00224'}; %have brain 

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspine leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

HFC=1;
fband=[10 35];

%load structure with positions
load('C:\Users\mspedden\Documents\brainspine_save\intersection_pos_M1.mat') %positions that intersect with anatomical cluster

% figure;
% ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.8 0.8], 'facealpha', 0.4, 'edgecolor', 'none');
% hold on;
% plot3(intersection_pos(:,1), intersection_pos(:,2), intersection_pos(:,3), ...
%       'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'LineWidth', 2);



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

    %remove bad channels
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

 brainlabs=[];
    for f = 1:length(grad_mm.label)
        if grad_mm.chanpos(f,2) > 200
            brainlabs = [brainlabs grad_mm.label(f)];
        end
    end

    %% new brain only dataset and grad

brainidx = find(grad_mm.chanpos(:,2) > 200);
braingrad = grad_mm;
braingrad.chanori = grad_mm.chanori(brainidx, :);
braingrad.chanpos = grad_mm.chanpos(brainidx, :);
braingrad.chantype = grad_mm.chantype(1:length(brainidx));
braingrad.chanunit = grad_mm.chanunit(1:length(brainidx));
braingrad.coilori = grad_mm.coilori(brainidx, :);
braingrad.coilpos = grad_mm.coilpos(brainidx, :);
braingrad.label = grad_mm.label(brainidx);
braingrad.tra = grad_mm.tra(brainidx, brainidx);

%figure; ft_plot_sens(braingrad)

braindat=ftdat;

mesh_brain.unit='mm';

%single shell volume conductor
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, mesh_brain);

% calculate forward
cfg                     = [];
cfg.sourcemodel         = sources_brain;
cfg.headmodel           = vol;
cfg.grad                = braingrad;
cfg.reducerank          = 'no';
LF = ft_prepare_leadfield(cfg);


%% lcmv

cfg                   = [];
cfg.covariance        = 'yes';
tlock = ft_timelockanalysis(cfg, braindat);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = vol;
cfg.sourcemodel.LF      = LF;
cfg.grid                = sources_brain;
cfg.unit                = LF.unit;
cfg.lcmv.keepfilter     = 'yes';
source_idx = ft_sourceanalysis(cfg, tlock);


% --- Map ROI positions to source grid ---
idx_roi = dsearchn(source_idx.pos, intersection_pos);
idx_roi = unique(idx_roi);

% keep only valid inside points
idx_roi = idx_roi(source_idx.inside(idx_roi));

roi_center = mean(source_idx.pos(idx_roi,:), 1);
d = sqrt(sum((intersection_pos - roi_center).^2, 2));
R = max(d);

d_all = sqrt(sum((source_idx.pos - roi_center).^2, 2));
roi_idx_radius = find(d_all <= R & source_idx.inside);

figure;
ft_plot_mesh(mesh_brain, ...
    'facecolor', [0.8 0.8 0.8], ...
    'facealpha', 0.7, ...
    'edgecolor', 'none');
hold on;

% plot ROI source points
plot3(source_idx.pos(roi_idx_radius,1), ...
      source_idx.pos(roi_idx_radius,2), ...
      source_idx.pos(roi_idx_radius,3), ...
      'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');



cfg = [];
cfg.pos          = roi_center;   
cfg.radius       = R;           
cfg.method       = 'svd';
cfg.numcomponent = 1;

VE = ft_virtualchannel(cfg, braindat, source_idx);


savename=sprintf('sub%s_VE_brain_forspectra',sub);
save(fullfile(save_dir,savename),"VE")

end %subjects




% beamformer=source_idx.avg.filter{1};
% 
% chansel  = tlock.label(1:end-1); % MEG sensors
% chanindx = 1:length(chansel);    %hard coded    
% 
% coh_data = [];
% coh_data.label = {'coh_x', 'coh_y', 'coh_z'};
% coh_data.time = braindat.time;
% 
% %this gives for each orientation
% for i=1:length(braindat.trial)
%   coh_data.trial{i} = beamformer * braindat.trial{i}(chanindx,:);
% end
% 
% timeseries = cat(2, coh_data.trial{:});
% [u1, s1, v1] = svd(timeseries, 'econ');
% 
% virtualchannelbrain = [];
% virtualchannelbrain.label = {'motor'};
% virtualchannelbrain.time = braindat.time;
% 
% for k = 1:length(braindat.trial)
%   virtualchannelbrain.trial{k}(1,:) = u1(:,1)' * beamformer * braindat.trial{k}(chanindx,:);
% end
% 
% savename=sprintf('brain_VC_%s',sub);
% save(fullfile(save_dir,savename), 'virtualchannelbrain')
% error('stop')
% 
% % cfg = [];
% % cfg.channel = 'EXG1';
% % emgdata = ft_selectdata(cfg, braindat);
% % 
% % cfg = [];
% % combineddata = ft_appenddata(cfg, virtualchannelbrain,virtualchanneldata);
% % 
% % cfg            = [];
% % cfg.output     = 'fourier';
% % cfg.method     = 'mtmfft';
% % cfg.foilim     = [5 45];
% % cfg.tapsmofrq  = 3;
% % cfg.keeptrials = 'yes';
% % freq = ft_freqanalysis(cfg, combineddata);
% % 
% % cfg = [];
% % cfg.method = 'coh';
% % coherence = ft_connectivityanalysis(cfg, freq);
% % figure; plot(coherence.freq,(squeeze(coherence.cohspctrm(1,:,:))'))
% % ylim([0 0.15])
% 

