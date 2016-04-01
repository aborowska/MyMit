n = 1000;
df = 1;
t_sample = random('t',df,[n,1]);
L100 = length(t_sample(abs(t_sample)<100));
hist(t_sample(abs(t_sample)<100), 50)

t_pdf = @(v,x) gamma((v+1)/2)/(sqrt(v*pi)*gamma(v/2)).*(1+(x.^2)/v).^(-(v+1)/2);
n_pdf = @(x) exp(-x.^2/2)./sqrt(2*pi);

x= -7:0.01:7;
t_v3 = t_pdf(3,x);
figure(2)
plot(x,t_v3)
title('Profit/loss density')
VaR_5 =  quantile(t_v3,0.05);
I = find(t_v3(t_v3<VaR_5));
max(I)
hold on
plot(x(1,50),0,'ko','MarkerFaceColor','r','MarkerSize',6);
plot(x(1,10),0,'ko','MarkerFaceColor','b','MarkerSize',6);
hold off
legend('Density','VaR','ES')


t_v100 = t_pdf(100,x);
t_v200 = t_pdf(200,x);


n_plot = n_pdf(x); 
hold on
plot(x,t_v1, 'b')
plot(x,t_v100,'r')
plot(x, n_plot, 'k')
hold off

plot(x,t_v200-n_plot,'r')

% Caychy
x = -10:0.01:10;
c = 1./(pi*(1+x.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%
v = [1,2,5,10,100];
g = length(v);

t_sample = @(df) random('t',df,[n,1]);
T = arrayfun(@(df)t_sample(v(df)), 1:g, 'un', 0);
T = cell2mat(T);


for ii = 1:g
    figure(ii)
    [A,B] = hist(T(:,ii),50);
    scatter(B,A)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
df = 8;
phi = randn(n,1);
Y_inv = random('gam',df/2,2/df,[n,1]);
Y = 1./Y_inv;
X = sqrt(Y).*phi;
T = random('t',df,[n,1]);

figure(1)
hist(X, 50)


figure(2)
hist(T, 50)

%%%%%%%%%%%%%%%%%%%%%%%%%%
df = 2;
prec = random('gam',df/2,df/2,[n,1]);
sigma = (1./prec).^(0.5);
mu = zeros(n,1);
n_sample = random('norm',mu,sigma);


%%%%%%%%%%%%%%%%%%%%%%%%%%
chi2_sample =  random('chi2',df,[n,1]);
gamma_sample = random('gam',df/2,2,[n,1]);

figure(1)
hist(chi2_sample, 50)

figure(2)
hist(gamma_sample, 50)