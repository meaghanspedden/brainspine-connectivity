                
%% setup minimum norm source reconstruction
%'OP00220', 'OP00221', 'OP00225', 'OP00226'};
%subs={'OP00224'}; %_002

sub='OP00225';

geomfile = 'D:\MSST001\generic_merged\geoms.mat';


datwithEMGmerged = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);


load(geomfile)
D=spm_eeg_load(datwithEMGmerged);
grad_mm=D.sensors('MEG');
ftdat=spm2fieldtrip(D);

channames=D.chanlabels;


    spinelabs={};

    for f = 1:length(grad_mm.label)
        if grad_mm.chanpos(f,2) < 200 
            spinelabs = [spinelabs grad_mm.label(f)];
        end

    end

    %% new sc only dataset and grad

[~,spineidx]=find(ismember(spinelabs, grad_mm.label)); %local idx grad struct
spinegrad=grad_mm;
spinegrad.chanori=grad_mm.chanori(spineidx,:);
spinegrad.chanpos=grad_mm.chanpos(spineidx,:);
spinegrad.chantype=grad_mm.chantype(1:length(spineidx));
spinegrad.chanunit=grad_mm.chanunit(1:length(spineidx));
spinegrad.coilori=grad_mm.coilori(spineidx,:);
spinegrad.coilpos=grad_mm.coilpos(spineidx,:);
spinegrad.label=grad_mm.label(spineidx);
spinegrad.tra=grad_mm.tra(spineidx,spineidx);


%only spinal cord channels
cfg=[];
cfg.channel=spinelabs;
spinedat=ft_selectdata(cfg,ftdat);

src_mm=sources_center_line;


                cfg                     = [];
                cfg.method              = 'infinite';
                cfg.siunits=1;
                cfg.grad=spinegrad; %mm
                cfg.conductivity = 1;

                vol                     = ft_prepare_headmodel(cfg,mesh_torso);
                vol.type                = 'infinite_currentdipole';
                vol.unit                = 'mm';

                % calculate forward
                cfg                     = [];
                cfg.sourcemodel         = src_mm;
                cfg.headmodel           = vol;
                cfg.grad                = grad_mm;
                cfg.reducerank          = 'no';
                %cfg.normalize           = 'yes';
                %cfg.normalizeparam     = 0.5;

                LF = ft_prepare_leadfield(cfg);

            cfg=[];
            cfg.output     = 'powandcsd';
            cfg.method     = 'mtmfft';
            cfg.foilim     = [10 35];
            cfg.tapsmofrq  = 2;
            cfg.keeptrials = 'no';
            freqdat=ft_freqanalysis(cfg,spinedat);

            %filter across both conditions
                cfg=[];
                cfg.grid = src_mm;
                cfg.headmodel=vol;
                cfg.sourcemodel.leadfield=LF;
                cfg.dics.keepfilter='yes';
                cfg.method = 'dics'; %depth normalised
                source_all = ft_sourceanalysis(cfg,freqdat); %gives one val pr source pt per freq and ori

               % figure; plot(src_mm.pos(:,2), mean(source_stat.avg.pow,2),'.-')
               

%% separate conditions

statidx=find(spinedat.trialinfo==1);
restidx=find(spinedat.trialinfo==2);


cfg=[];
cfg.trials=statidx;
statdat=ft_selectdata(cfg,spinedat);

cfg.trials=restidx;
restdat=ft_selectdata(cfg,spinedat);

            cfg=[];
            cfg.output     = 'powandcsd';
            cfg.method     = 'mtmfft';
            cfg.foilim     = [10 35];
            cfg.tapsmofrq  = 2;
            cfg.keeptrials = 'no';
            freqstat=ft_freqanalysis(cfg,statdat);
            freqrest=ft_freqanalysis(cfg,restdat);

               %now apply filter to each condition
                cfg=[];
                cfg.grid = src_mm;
                cfg.headmodel=vol;
                cfg.sourcemodel.leadfield=LF;
                cfg.dics.filter=source_all.avg.filter;
                cfg.method = 'dics'; %depth normalised
                source_stat = ft_sourceanalysis(cfg,freqstat);

                source_rest = ft_sourceanalysis(cfg,freqrest);

                figure; plot(src_mm.pos(:,2), mean(source_stat.avg.pow,2),'.-')
                hold on; plot(src_mm.pos(:,2), mean(source_rest.avg.pow,2),'.-')
                legend({'static', 'rest'})

