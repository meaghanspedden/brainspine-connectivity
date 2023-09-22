

% parameters

n = 3000; % rows

df1=31;  % columns matrix 1

df2=31;  % columns matrix 2

 

% data

a= randn(n,df1);

b= randn(n,df2);

 

% cca simulated

Ein = orth(a);

Eout = orth(b);

[u,ss,v] = svd(Ein'*Eout);

 

% cca predicated

c1=df1/(n);

c2=df2/(n);

r1= sqrt(c1*c2/((1-c1)*(1-c2)));

sst =r1*(1-c1+c1/r1)*(1-c2+c2/r1);

 

% simulated and predicted, squared canonical correlations

[ss(1,1)^2 sst]