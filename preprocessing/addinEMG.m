
function [dD,EMGsamples,EMGdata]=addinEMG(D,dD,fullEMGname,EMGtriglatency,OPMEMGsync,OPMchansreplace,EMGchannames);
%function dD=addinEMG(D,dD,EMGstruct,EMGtriglatency,OPMEMGsync,OPMchanreplace,EMGchanname);



if size(OPMchansreplace,1)~=size(EMGchannames,1),
    error('OPM and EMG labels not the same size')
end;


synchTrigIdx = find(strcmp(D.chanlabels,OPMEMGsync)); 
%figure; plot(rawData.trial{1}(synchTrigIdx,:),'k')

%find when trig comes
eventCont           = D(synchTrigIdx,:,1);
time                = D.time;
sample              = 1:length(time);
% convert to binary
eventDisc           = (eventCont > 3);

% find first sample. first find differences
tmp                 = eventDisc - [0,eventDisc(1:end-1)];
tmp2                = (tmp == 1);
events              = sample(tmp2)';
events              = round(events);

%because it is downsampled
event = round((events/(D.fsample/dD.fsample)));

%trigsample=event; %% min(find(dtrigchan>max(dtrigchan)/2))
%EMGstartsample=trigsample

EMGstartsample=event; %

for f=1:size(OPMchansreplace,1),
    [~, EMGstruct]= readSpikeData(fullEMGname,deblank(EMGchannames(f,:)));
    %EMGdata(f,:)=detrend(EMGstruct.trial{1});
    d1=EMGstruct.trial{1};
    EMGdata(f,:)=d1(EMGtriglatency*EMGstruct.fsample:end); %% discard pre-trig samples
    EMGdata(f,:)=detrend(EMGdata(f,:));
end;

if dD.fsample~=EMGstruct.fsample,
    error('downsampled datasets must have same sample rate');
end;


nEMGsamples=length(EMGdata);

EMGsamples=EMGstartsample:EMGstartsample+nEMGsamples-1;
if max(EMGsamples)>dD.nsamples,
    EMGsamples=EMGstartsample:dD.nsamples;
    nEMGsamples=length(EMGsamples);
end;


for rp=1:size(OPMchansreplace,1),
    OPMchanreplace=deblank(OPMchansreplace(rp,:));
    replaceind=dD.indchannel(OPMchanreplace);
 
    if isempty(replaceind)
        error('could not find channel to replace');
    end;
    dD(replaceind,:,1)=dD(replaceind,:,1).*0;
    dD(replaceind,EMGsamples,1)=EMGdata(rp,1:nEMGsamples);
    dD = chanlabels(dD,replaceind,deblank(EMGchannames(rp,:)));
    dD = units(dD,replaceind,'uV');
    dD = chantype(dD,replaceind,'EMG');

end; %for replaceind
dD.save;
PLOTON=0;
if PLOTON,
    figure;
    
    plot(dD.time,squeeze(dD(replaceind,:,1)),dD.time,dD(synchTrigIdx,:,1)*std(dD(replaceind,:,1)))
end;


