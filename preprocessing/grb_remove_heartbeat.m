function D=grb_remove_heartbeat(D,heartep,chind,megind,BALANCE);
%function D=grb_remove_heartbeat(D,chind,BALANCE);
%% use  back channels to get estimate of beat timing
% if BALANCE THEN CORRECT montage in file


remind=setdiff(megind,chind); %% indices of good chans not on lower back

cdata_lb=heartep(chind,:)*heartep(chind,:)';
cdata_rem=heartep(remind,:)*heartep(remind,:)';
[u_lb,s_lb,v]=svd(cdata_lb);
[u_rem,s_rem,v]=svd(cdata_rem);
figure;
dsum_lb=cumsum(diag(s_lb))./sum(diag(s_lb));
dsum_rem=cumsum(diag(s_rem))./sum(diag(s_rem));
plot(1:length(diag(s_lb)),dsum_lb,1:length(diag(s_rem)),dsum_rem)
legend('lb','rem')
Ncomp_lb=min(find(dsum_lb>0.99));
Ncomp_rem=min(find(dsum_rem>0.99));
fprintf('\n Removing %d heart components from data and %d from remainder chans',Ncomp_lb,Ncomp_rem);

%% project heart out of back chans
S=[];
S.D=D;
S.X=u_lb(:,1:Ncomp_lb);
S.channels=chind;
S.balance=BALANCE; %% correct lead fields
D = hfo_project_out_comps(S);

%% project heart out of remainder channels
S=[];
S.D=D;
S.X=u_rem(:,1:Ncomp_rem);
S.channels=remind;
S.balance=0; %% do not correct lead fields (as no grad structure defined)
D = hfo_project_out_comps(S);
D.save;



