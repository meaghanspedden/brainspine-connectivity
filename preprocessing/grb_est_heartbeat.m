function [heartep,beatlen,megind]=grb_est_heartbeat(D,chind,hbcomp,beatlen,megind);


%% use  svd to get estimate of heartbeat
if nargin<4
    beatlen=[];
end

if nargin<5
    megind=[];
end

if isempty(megind)
    megind=setdiff(D.indchantype('MEG'), D.badchannels);
end

remind=setdiff(megind,chind); %% indices of good chans not on back

data=squeeze(D(chind,:,1));

cdata=data*data';

[u,s,v]=svd(cdata);

heartest=abs(data'*u(:,hbcomp));

thresh=std(heartest)*4;

figure
plot(D.time,abs(heartest),D.time,thresh*ones(size(D.time)));

%beattrigs=[0 ;diff(heartest>thresh)];
% if isempty(beatlen)
%     beatlen=2*floor((median(diff(find(beattrigs==1))))/2); % even
% end
%% now make average of all MEG channels based around beat

[pks,  htrigs]=findpeaks(heartest,'MinPeakDistance',400,'MinPeakHeight',thresh);
findpeaks(heartest,'MinPeakDistance',400,'MinPeakHeight',thresh);
hold on

pct95pks=prctile(pks,95);
outliers=find(pks>pct95pks);
vline(htrigs(outliers))
htrigs(outliers)=[];

%heart beat length is median diff between peaks
if isempty(beatlen)
    beatlen=2*floor((median(diff(htrigs)))/2); % even
end

hbfreq=length(htrigs)/D.time(end)*60;
fprintf('est heart rate %g bpm\n',hbfreq)

heartep=[];

for f=1:length(htrigs)

    if (((htrigs(f)+beatlen)<length(heartest)) && ((htrigs(f)-beatlen)>0)) 
        heartep(:,:,f)=D(megind,htrigs(f)-beatlen/2:htrigs(f)+beatlen/2-1,1);
    end

end

heartep=mean(heartep,3);


f=figure;
plot(heartep')
legend('estimate of heartbeat (over all chans)')
waitfor(f)

