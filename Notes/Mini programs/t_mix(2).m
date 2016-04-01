n = 6000;
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