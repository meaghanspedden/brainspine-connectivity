%% SET UP FORWARD MODEL
N=1000; %% time points
K=1; %% latent factors
M=10; %% measurement channels
A=randn(M,K); %% forward model
e=1; %% noise level
s=zeros(K,M); % latent factor
for k=1:K,
    s(:,1:N)=ones(1,N); %sin((1:N)/20); %% time course of latent factor
    s(:,1:N/2)=0;
end;
x=zeros(N,M); %% data
x=A*s+randn(M,K).*e; %% eqn 1 from Haufe paper
x=x-repmat(mean(x')',1,size(x,2)); %% dc correct
plot(1:N,s,'k',1:N,x,'r');
legend('latent time course','measured time course');

%% NOW TRY AND INTERPRET
XD=ones(N,1); %% design matrix
XD(N/2:end)=-XD(N/2:end);
Yred=x'; %% measured data
[CVA] = spm_cva(Yred,XD);


covdat=cov(Yred);
Aest=covdat*CVA.V*inv(cov(CVA.v)); %% from george+ Meaghan also eqn 8 from Haufe paper
uAest=Aest./sqrt(dot(Aest,Aest)); %% remove scaling factor
uA=A./sqrt(dot(A,A));
if corr(uA,uAest)<0, %% can get a sign flip
    uAest=-uAest;
end; 
plot(1:M,uA,'-gx',1:M,uAest,'ro');
xlabel('channel')
ylabel('estimate of weight')
legend('red=est, green=true')


