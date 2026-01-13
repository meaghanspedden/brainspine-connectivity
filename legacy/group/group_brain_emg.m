% summarize brain-EMG beta power modulations across participants

clear all
close all

subs = {'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219', ...
         'OP00221', 'OP00224', 'OP00225', 'OP00226'};

all_betacoh = [];
all_powdiff = [];

load('D:\brainspineconnectivity\topoplot\brain_layout.mat')

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    fname = fullfile(save_dir, ['brainEMGmod_subject_' sub '.mat']);
    if exist(fname, 'file')
        dat = load(fname);

        if isempty(all_powdiff)
            n_chans=length(dat.log_power_diff);
            all_powdiff = nan(length(subs), n_chans);
            all_betacoh = nan(length(subs), n_chans-1);

        end

        all_powdiff(i,:) = dat.log_power_diff;
        all_betacoh(i,:) = dat.betacoh;

    else
        warning('File not found for %s', sub);
    end
end

%% Group-level t-stats for brain beta modulation
[~, p, ~, stats] = ttest(all_powdiff);
group_braint = stats.tstat(1:end-1);


figure
topobrain(brainlay, group_braint')
title('Group-level t-stat: Brain beta')
hold on
scatter(brainlay.pos(3,1),brainlay.pos(3,2), 'ro')


meancoh=mean(all_betacoh,1);

figure
topobrain(brainlay, meancoh')
title('Mean beta pow corr')
hold on
scatter(brainlay.pos(3,1),brainlay.pos(3,2), 'ro')
