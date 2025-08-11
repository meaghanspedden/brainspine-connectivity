%% Beamformer
%outputs time series for each of the 3 orientations and each sourcepoint
%create_geoms (1. setup_source; 2. create_leadfields 3. beamforming)
clear all
close all

sub='OP00215';

addpath('D:\msg_source_recon\')

save_dir = fullfile('D:\MSST001', [sub '_contrast']);
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
subject_dir = fullfile('D:\MSST001', [sub '_static']); %rest and static will be same here

filename_rest = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['prestoe1000mspddfflo45hi45hfcrest_001_array1.mat']);

filename_static = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pstaticoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);


%% fixed across participants
brainchan_labels={'H1', 'H6', 'B3', 'H8', 'G3', 'B6', 'H7', 'H4', 'G4', 'B4', 'B2', 'H2', 'B7', 'B8', 'B5', 'B1'};

brainlabs=labelToIndex(brainchan_labels);

%%

% SETUP AND PARAMS
covtime = [-Inf Inf];  % Covariance time window
hstep = 10;  % Grid spacing for mesh sampling brain

geomfile = fullfile(subject_dir, 'geoms.mat');
castforward = load(geomfile);

% PLOT MESH & SOURCE POSITIONS
mtorso = castforward.mesh_torso;
spos = castforward.sources_center_line.pos;
Ncordpoints = length(spos);



%% BRAIN------------------------------------------------------------------
% Identify head region
headind = find(mtorso.vertices(:,2) > 250); %spos(1,2));
%plot3(mtorso.vertices(headind,1),mtorso.vertices(headind,2),mtorso.vertices(headind,3),'mo');

% Generate head mesh and sample inside
hvert = mtorso.vertices(headind,:);
dT = delaunayTriangulation(hvert);
hfaces = convexHull(dT);

headmesh=[];
headmesh.faces=hfaces;
headmesh.vertices=hvert;
headmesh.unit='mm';

M.faces = hfaces;
M.vertices = hvert;
M = gifti(M);
allhp = [];
isinside=[];


for xh = min(hvert(:,1)):hstep:max(hvert(:,1))
    for yh = min(hvert(:,2)):hstep:max(hvert(:,2))
        for zh = min(hvert(:,3)):hstep:max(hvert(:,3))
            allhp = [allhp; xh yh zh];

            if spm_mesh_inside(M, [xh yh zh])
                isinside=[isinside 1];
            else
                isinside=[isinside 0];
            end
        end
    end
end
isinside=logical(isinside);
plot3(allhp(isinside,1), allhp(isinside,2), allhp(isinside,3), '.') %source points in head area

headmesh2=[]; %fieldtrip format
headmesh2.pos=headmesh.vertices;
headmesh2.tri=headmesh.faces;
headmesh2.unit=headmesh.unit;

cfg = [];
cfg.method = 'headshape';
cfg.headshape = headmesh2;  
cfg.headshape.check = 'yes';  
mesh_checked = ft_prepare_mesh(cfg);

%brain volume conductor
cfg=[];
cfg.method='singleshell';
vol=ft_prepare_headmodel(cfg,mesh_checked);

%brain source
brainsrc=[];
brainsrc.pos=allhp;
brainsrc.inside=isinside';
brainsrc.unit='mm';

% LOAD  DATA--------------------------------------------------------
Drest = spm_eeg_load(filename_rest);
Dstatic = spm_eeg_load(filename_static);

%convert to FT for beamforming
ftdat_rest=spm2fieldtrip(Drest);
ftdat_static=spm2fieldtrip(Dstatic);

%append conditions to make common filter
alldat=ft_appenddata([],ftdat_static,ftdat_rest); %

%prep chans and idx
grad=Drest.sensors('MEG');
labs=Drest.chanlabels;
%can do this with ismemeber instead
brainchans = labs( ...
    contains(labs, brainlabs) & ...
    (contains(labs, 'X') | contains(labs, 'Y') | contains(labs, 'Z')));
brainchanidx=find(ismember(labs,brainchans));
brainchanidx = setdiff(brainchanidx, Drest.badchannels);

%get leadfields
cfg = [];
cfg.grad = grad;
cfg.headmodel = vol;
cfg.sourcemodel = brainsrc;
cfg.reducerank = 'no';
cfg.channel = labs(brainchanidx);
cfg.normalize = 'no';
leadfield_brain = ft_prepare_leadfield(cfg);

%insideindx=find(isinside);

figure; hold on
ft_plot_sens(grad, 'style', '*b');
ft_plot_headmodel(vol, 'facecolor', 'm','edgecolor','none'); alpha 0.2;
plot3(brainsrc.pos(:,1), brainsrc.pos(:,2),brainsrc.pos(:,3),'ko');

%get power for dics
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.keeptrials='yes';
cfg.tapsmofrq = 4;
cfg.foilim    = [15 30];
freqdat = ft_freqanalysis(cfg, alldat);
fdatstat=ft_freqanalysis(cfg,ftdat_static);
fdatrest=ft_freqanalysis(cfg, ftdat_rest);

%get filter using both conditions
cfg                  = [];
cfg.method           = 'dics';
cfg.sourcemodel      = leadfield_brain;
cfg.headmodel        = vol;
cfg.channel          = labs(brainchanidx);
cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.dics.keepfilter  ='yes';
source_brain         = ft_sourceanalysis(cfg, freqdat);


% %apply filter for single trials
cfg                  = [];
cfg.sourcemodel      = leadfield_brain;
cfg.sourcemodel.filter=source_brain.avg.filter;
cfg.headmodel        = vol;
cfg.channel          = labs(brainchanidx);
cfg.rawtrial     ='yes';
source_stat  = ft_sourceanalysis(cfg, fdatstat);
source_rest  = ft_sourceanalysis(cfg, fdatrest);

ntrials=min([length(source_stat.trial), length(source_rest.trial)]);

restmat=[];
statmat=[];

for k=1:ntrials

restmat(k,:)=source_stat.trial(k).pow(brainsrc.inside);
statmat(k,:)=source_rest.trial(k).pow(brainsrc.inside);

end

[H, P, CI, STATS] = ttest2(log(statmat), log(restmat)); % 1 t stat per source point

func=STATS.tstat;

cmin = min(func(:));
cmax = max(func(:));
clim = max(abs(func(:)));

[pkpos, peakind_pos] = max(func);  %indexed thru inside
[pkneg, peakind_neg] = min(func); 

if abs(pkneg) > pkpos
    fprintf('peak t is neg/n')
    idx2plot=peakind_neg;
else
    fprintf('peak t is pos/n')
    idx2plot=peakind_pos;
end

inside=brainsrc.pos(find(brainsrc.inside),:);
fig=figure; hold on;

ft_plot_mesh(mtorso,'FaceColor', [0.8 0.8 1.0], 'facealpha', 0.2, 'EdgeColor','none',...
    'AmbientStrength', 0.15); hold on
ft_plot_cloud(inside,func,'scalerad','no', 'cloudtype','surf','radius',3.5,'clim', [-clim clim])
xlabel('X'); ylabel('Y')
axis image
plot3(inside(idx2plot,1), inside(idx2plot,2),inside(idx2plot,3),'ro','MarkerSize',14)
colorbar

plot3(spos(:,1),spos(:,2),spos(:,3),'ro'); axis equal;

waitfor(fig)

 % Save results
 bffilename = fullfile(save_dir, ['bfdata_brain_', sub]);

 save(bffilename, 'statmat', 'brainsrc', 'func','idx2plot','mtorso');















% %make into 3d volume for nifti printing
% nX = length(unique(allhp(:,1)));
% nY = length(unique(allhp(:,2)));
% nZ = length(unique(allhp(:,3)));
%
% %map grid positions to volume indices
% xvals = unique(allhp(:,1));
% yvals = unique(allhp(:,2));
% zvals = unique(allhp(:,3));
%
% voxsize = [xvals(2)-xvals(1), yvals(2)-yvals(1), zvals(2)-zvals(1)];
% origin  = [xvals(1), yvals(1), zvals(1)];
%
% affine = [voxsize(1), 0,          0,          origin(1);
%          0,          voxsize(2), 0,          origin(2);
%          0,          0,          voxsize(3), origin(3);
%          0,          0,          0,          1];
%
% V = struct();
% V.dim = [nX nY nZ];
% V.dt = [16 0];  % float32
% V.mat = affine;
% V.descrip = 'Functional source data';
%
% for im=1:5 %length(source_brain_tr.trial)
% % Step 3: For each inside point, get indices in volume
% vol3D = nan(nX, nY, nZ);
%
% for i = 1:size(allhp,1)
%     if isinside(i)
%         xidx = find(xvals == allhp(i,1));
%         yidx = find(yvals == allhp(i,2));
%         zidx = find(zvals == allhp(i,3));
%         vol3D(xidx, yidx, zidx) = source_brain_tr.trial(im).pow(i);
%     end
% end
%
%
% V.fname = sprintf('trial_%g_%s.nii',im,analysis);
% cd(save_dir)
% spm_write_vol(V, vol3D);
%
%
% end



% tic
% brain_source=[]; %loop through src pts and trials
% for s=1:length(brainsrc.pos)
%
%     for t=1:length(ftdat.trial)
%         bf_filt = source.avg.filter{s};
%         brain_source(s,:,t)= bf_filt * ftdat.trial{t}(brainchanidx,:);
%     end
%
% end
% toc





   


