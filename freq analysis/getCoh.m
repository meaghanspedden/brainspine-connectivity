%% sensor level coherence


sub = 'OP00225';
concat = 1;  % Set to 0 to load only one run
run2use='001';

if concat %note that this has hb removed and has rest!!!

    filename = fullfile('D:\MSST001', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

    DAll=spm_eeg_load(filename);

else

    filename = fullfile('D:\MSST001', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        ['oe1000mspddfflo45hi45hfcstatic_' run2use '_array1.mat']);

    DAll=spm_eeg_load(filename);
end

%%
load('D:\brainspineconnectivity\brain_layout.mat')

Zbrainchans=brainlay.label;

%brainchans=labelToIndex(brainchan_labels);
% todel=find(startsWith(brainchans, {'39', '29', '19', '51', '58', '59', '49'}));
% brainchans(todel)=[];
%Zbrainchans = cellfun(@(x) ['Z' x], brainchans, 'UniformOutput', false);

emgidx=find(contains(DAll.chanlabels,'EXG1'));

ftdat=spm2fieldtrip(DAll);
if concat %only static condition
    cfg=[];
    cfg.trials=find(ftdat.trialinfo==1);
    ftdat=ft_selectdata(cfg,ftdat);
end


cfg=[];
cfg.channel=Zbrainchans;
opmdat=ft_selectdata(cfg, ftdat);

cfg=[];
cfg.channel='EXG1';
cfg.rectify ='yes';
cfg.detrend='yes';
EMG=ft_preprocessing(cfg, ftdat);

alldat=ft_appenddata([], opmdat, EMG); %use this for Neurospec


cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 45];
cfg.tapsmofrq  = '2';
cfg.keeptrials = 'yes';
% cfg.pad='maxperlen';

freq_dat   = ft_freqanalysis(cfg, alldat);



seed='EXG1';
brain=Zbrainchans;
tmpcmb=repmat({seed}, length(brain),1);
combs=[brain tmpcmb];

cfg            = [];
cfg.method     = 'coh';
cfg.complex    ='absimag';
cfg.channelcmb = combs;
coh      = ft_connectivityanalysis(cfg, freq_dat);

%     fig=figure;
%     cfg=[];
%     ft_connectivityplot(cfg,coh)
%     waitfor(fig)

% for k=1:length(brain)
%     f=figure;
%     plot(coh.freq, coh.cohspctrm(k,:),'LineWidth',2)
%     title(sprintf('%s',Zchans{k}))
%     waitfor(f)
% end


figure; plot(coh.freq, coh.cohspctrm)

betafreqs=find(coh.freq >=15 & coh.freq<=30);
betacoh=mean(coh.cohspctrm(:,betafreqs),2);

% braingrad=grad;
% braingradidx=find(contains(grad.label,Zbrainchans));
% 
% braingrad.chanori=grad.chanori(braingradidx,:);
% braingrad.chanpos=grad.chanpos(braingradidx,:);
% braingrad.coilori=grad.coilori(braingradidx,:);
% braingrad.coilpos=grad.coilpos(braingradidx,:);
% braingrad.label=grad.label(braingradidx);
% braingrad.chantype=grad.chantype(braingradidx);
% braingrad.chanunit=grad.chanunit(braingradidx);
% braingrad.tra=grad.tra(braingradidx,braingradidx);

% figure; ft_plot_mesh(geoms.mesh_torso, 'facealpha',0.7, 'edgealpha', 0.2)
% hold on
% ft_plot_sens(grad,'label','on')
% view(179.2635,  -6.3792)

figure
subplot(211)
topobrain(brainlay, betacoh)

% Sum coherence across frequencies for each channel
coh_sum = sum(coh.cohspctrm, 2);  % sum over freqs (dim 2)

% Find (local) index of maximum sum
[~, max_idx] = max(coh_sum);
brainlay.label(max_idx)

subplot(212)
plot(coh.freq, coh.cohspctrm(max_idx,:))

% brainchan2use=labels(max_idx);
% fprintf('brainchan max coh:%s/n', brainchan2use)

%neurospec coh

% addpath('D:\neurospec20')
% addpath(genpath('D:\neurospec211NEW'))
%
% contdat=[];
% for k=1:length(alldat.trial)
%     contdat=[contdat alldat.trial{k}];
% end
%
%
% opt_str='t2 M5';
% seg_pwr=10;
%
% cohmat=[];
% f11=[];
% f22=[];
%
% EMG=contdat(emgidx,:);
%
% for cc=1:length(brainchanidx)
%     dat1=contdat(brainchanidx(cc),:)';
%     [f,t,cl,sc] = sp2a2_R2_mt(dat1,EMG',alldat.fsample,seg_pwr,opt_str);
%     cohmat(:,cc)=f(:,11)+f(:,12);
%
%     %f(:,4);
%     f11(:,cc) = f(:,2);
%     f22(:,cc) = f(:,3);
% end
%
% % for k=1:length(brainchanidx)
% %     fig=figure;
% %     plot(f(:,1), cohmat(:,k))
% %     xlim([0 100])
% %     waitfor(fig)
% % end
%
% medcoh=median(cohmat,2);
%
% figure; plot(f(:,1), cohmat)
% hold on
% plot(f(:,1), medcoh, 'k', 'LineWidth',2)
% xlim([0 100])