
%% source recon experiments in fieldtrip using BEM leadfields spine

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

%sub='OP00212';
subs = {'OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_realistic_cervical_brainspine_fine.mat');

HFC=1;
LFop='spine';
rectify=1;
fband=[10 35];

subjResults=struct();

for ss=1:length(subs)

datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');


load(geomfile)
D=spm_eeg_load(datwithEMGmerged);
grad_mm=D.sensors('MEG');
ftdat = spm2fieldtrip(D);

badchans=D.chanlabels(D.badchannels);

%remove bad channels here.
cfg=[];
cfg.channel=setdiff(ftdat.label,badchans);
ftdat=ft_selectdata(cfg,ftdat);

%%
if rectify
    cfg=[];
    cfg.rectify='yes';
    cfg.channel='EXG1';
    ftdatr=ft_preprocessing(cfg,ftdat);

    for k=1:length(ftdat.trial)
        ftdat.trial{k}(end,:)=ftdatr.trial{k}; %ftdat has rectified emg
    end

end


[Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'complex_bone_fields_cervical_brainspine_finegrid'), LFop);

nsourcepoints = size(Gx,1);
nchannels     = size(Gx,2);
spchanidx=find(sources_cent.pos(:,2) < 200);

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
cfg.foilim     = fband;
cfg.tapsmofrq  = 1;
cfg.keeptrials='yes';
freqdat_tr=ft_freqanalysis(cfg,ftdat); %trial wise freq

cfg=[];
cfg.avgoverfreq='yes';
freqdat_tr=ft_selectdata(cfg,freqdat_tr);%need to average over frequency as input to permutation


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

[nTrials,~] = min([length(statidx) length(restidx)]);

cfg=[];
cfg.trials=statidx(1:nTrials);
statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

cfg.trials=restidx(1:nTrials);
restdat=ft_selectdata(cfg,freqdat_tr);



%% get struct with mean freq data for constructing bf filter
cfg=[];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = fband;
cfg.tapsmofrq  = 1;
cfg.keeptrials = 'no';
freqdat=ft_freqanalysis(cfg,ftdat);

cfg=[];
cfg.avgoverfreq='yes';
freqdat=ft_selectdata(cfg,freqdat);

%% DICS coherence filter constructed on both conditions

cfg=[];
cfg.grid = sources_cent;
cfg.headmodel=dummyvol;
cfg.sourcemodel.leadfield=Lf;
cfg.dics.keepfilter='yes';
cfg.dics.lambda=10;
cfg.method = 'dics'; 
cfg.refchan='EXG1';

coh_source = ft_sourceanalysis(cfg,freqdat); 


%% apply to each condition separately
cfg=[];
cfg.grid = sources_cent;
cfg.headmodel=dummyvol;
cfg.sourcemodel.leadfield=Lf;
cfg.dics.filter=coh_source.avg.filter;
cfg.dics.lambda=10;
cfg.method = 'dics';
cfg.refchan='EXG1';

% source_stat = ft_sourceanalysis(cfg,statdat); %struct per condition with single trials
% source_rest = ft_sourceanalysis(cfg,restdat);

cfg.permutation = 'yes';
cfg.numpermutation=500;
source_perm = ft_sourceanalysis(cfg, statdat, restdat);


nPerm = numel(source_perm.trialA);
cohDiff_perm = zeros(nsourcepoints, nPerm);

for i = 1:nPerm
    cohA_perm = source_perm.trialA(i).coh;
    cohB_perm = source_perm.trialB(i).coh;
    cohDiff_perm(:, i) = cohA_perm - cohB_perm;
end

thr95 = prctile(cohDiff_perm, 95, 2);

%% DICS results
figure('Color','w'); hold on
diffCoh = source_perm.avgA.coh - source_perm.avgB.coh;
hLine = plot(sources_cent.pos(:,2), diffCoh, 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
hThresh = plot(sources_cent.pos(:,2), thr95, 'k--', 'LineWidth', 1.5);
aboveIdx = diffCoh > thr95;
xFill = [sources_cent.pos(aboveIdx,2); flipud(sources_cent.pos(aboveIdx,2))];
yFill = [thr95(aboveIdx); flipud(diffCoh(aboveIdx))];
hFill = fill(xFill, yFill, [0.9 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('Cranial–caudal distance (mm)', 'FontWeight', 'bold', 'FontSize', 12)
ylabel('Contraction − Rest', 'FontWeight', 'bold', 'FontSize', 12)
grid on; box on
set(gca, 'FontSize', 12)
legend([hLine hThresh hFill], {'Difference', '95th percentile', 'Significant'}, 'Location', 'best')
%% source structure with coh difference

coh_diff=source_perm.avgA.coh-source_perm.avgB.coh;
source_diff=coh_source;
source_diff.avg.coh=coh_diff;

cfg = [];
cfg.parameter = {'coh'};
spine_int=ft_sourceinterpolate(cfg,source_diff, mesh_wm);


%% add a mask
source_mask=coh_source; %copy
source_mask.avg.pow = coh_diff > thr95;

%interpolate the mask
cfg = [];
cfg.parameter = 'pow';          
cfg.interpmethod = 'nearest';   
source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_wm);
spine_int.mask=source_mask_int.pow;


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
figure
cfg = [];
cfg.figure='gcf';
cfg.method = 'surface';
cfg.funparameter = 'coh';
cfg.funcolormap = 'magma';
cfg.funcolorlim='zeromax';
cfg.projmethod = 'nearest';
cfg.surffile = mesh_wm;   
ft_sourceplot(cfg, spine_int); 
hold on
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.2, 'edgecolor', 'none');
ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.07, 'edgecolor', 'none');
ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.1, 'edgecolor', 'none');
view( -250, -1)
camlight

%% with mask
figure
cfg = [];
cfg.figure='gcf';
cfg.method = 'surface';
cfg.funparameter = 'coh';
cfg.funcolormap = 'magma';
cfg.maskparameter = 'mask';
cfg.maskstyle='opacity';
%cfg.funcolorlim='minzero';
%cfg.funcolorlim = [-.034 -.03];
cfg.projmethod = 'nearest';
cfg.surffile = mesh_wm;   
ft_sourceplot(cfg, spine_int); 
hold on
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.2, 'edgecolor', 'none');
ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.07, 'edgecolor', 'none');
ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.1, 'edgecolor', 'none');
view( -250, -1)
camlight


subjResults(ss).coh_diff=coh_diff;        % source_perm.avgA.coh - avgB.coh
subjResults(ss).thr95=thr95; 

end

nSubjects = length(subs);
sig_pos = false(nSubjects,1);

for ss = 1:nSubjects
    diffCoh = subjResults(ss).coh_diff;   % source_perm.avgA.coh - avgB.coh
    thr95    = subjResults(ss).thr95_pos; % 95th percentile from permutation

    if any(diffCoh > thr95)
        sig_pos(ss) = true;
    end
end

fprintf('Permutation: %d/%d subjects show a positive effect above threshold\n', ...
    sum(sig_pos), nSubjects);


%%
% figure;
% 
% % Brain mesh: semitransparent red/pink
% ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.5, 'edgecolor', 'none');
% hold on;
% 
% % Torso mesh: transparent blue shell
% ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');
% 
% % White matter mesh: solid white/grey core
% ft_plot_mesh(mesh_wm, 'facecolor', [0.9 0.9 0.9], 'facealpha', 1, 'edgecolor', 'none');
% 
% % Bone mesh: light beige, semitransparent
% ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
% 
% axis equal; 
% camlight; lighting gouraud
% %%
% figure
% ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');
% hold on
% plot3(sources_combined.pos(:,1), sources_combined.pos(:,2), sources_combined.pos(:,3),'b.')
% 
