function D=remove_heartbeat(D,heartep,megind,BALANCE)
%function D=grb_remove_heartbeat(D,chind,BALANCE);

% if BALANCE THEN CORRECT montage in file


cdata=heartep*heartep';

[u,s,v]=svd(cdata);

figure;
dsum=cumsum(diag(s))./sum(diag(s));
figure
plot(1:length(diag(s)),dsum)

Ncomp=min(find(dsum>0.99));
fprintf('\n Removing %d heart components from data',Ncomp);

%% project heart out of back chans
S=[];
S.D=D;
S.X=u(:,1:Ncomp);
S.channels=megind;
S.balance=BALANCE; %% correct lead fields
D = hfo_project_out_comps(S);

D.save;



