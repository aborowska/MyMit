ata = csvread('GSPC_ret.csv')
data = 100*data;
plot(data)


% QERMit ARCH 
I_arch = find(data<=-5.5);
data_arch = data(1:576,1);
S = var(data_arch);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mit=mit_init;
mu=mit.mu;
d=length(mu);
Sigma=reshape(mit.Sigma,d,d);

t1 = theta_new(:,1);
t2 = theta_new(:,2);


t1 = theta(:,1);
t2 = theta(:,2);

theta = csvread('theta.csv');
w = csvread('w.csv');

lnkR = csvread('lnk.csv');
lndR = csvread('lnd.csv');

n=1000;
figure(1)
df=1;
% draw = rmvt(mu,Sigma,df,n);
draw = rmvt([7,5],[1,0.5;0.5,1],4,n);

t1= draw(:,1);
t2= draw(:,2);
scatter(t1,t2)

figure(2)
df=100;
draw2 = rmvt(mu,Sigma,df,n);
t21= draw2(:,1);
t22= draw2(:,2);
scatter(t21,t22)

figure(3)
draw3=mvnrnd(mu,Sigma,n);
t31= draw3(:,1);
t32= draw3(:,2);
scatter(t31,t32)

C = cov(t1,t2);
rho = C(1,2)/sqrt(C(1,1)*C(2,2));


mu = [7, 5];
Sigma = [1, 1/2; 1/2, 1]; 
df = 4;
n = 4e3;

tic
figure(1)
draw= rmvt(mu,Sigma,df,n);
t1= draw(:,1);
t2= draw(:,2);
scatter(t1,t2)
time = toc;

figure(2)
draw2= rmvt_new(mu,Sigma,df,n);
t21= draw2(:,1);
t22= draw2(:,2);
scatter(t21,t22)


cov(draw)
cov(draw2)
mean(draw)
mean(draw2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = repmat([1,4,7,2,5,8,3,6,9],10,1)
AA = A * kron(eye(3),ones(3,1)); % the sum in the denominator of w(theta)
BB = A * repmat(diag(ones(3,1)), 3, 1)


theta_R = csvread('theta_R.csv');
lnd_R = csvread('lnd_R.csv');








sp=[(2514:-1:1)',AdjClose]
[B,I]=sort(sp(:,1))
spp=sp(I,:)
spp=spp(:,2)
I=1:2513;
sp_y=100*log(spp(I+1)./spp(I));