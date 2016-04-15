% addpath(genpath('include/'));

clear all;
close all;
%% Sampling of all the candidates at once
p = @(xx) 6*(xx-xx.^2);
q = @(xx) 1;
N = 15000;
BurnIn = 5000;

x_cand = rand(N,1);

u = rand(N,1);
x_chain = zeros(N,1);
x_chain(1,1) = x_cand(1,1);

acc = @(x0,x1) min(1,exp(log(p(x1)) - log(p(x0)) + log(q(x1)) - log(q(x0))));
p_val = p(x_cand);
q_val = q(x_cand);
lnw_v = log(p_val) + log(q_val);
acc2 = @(x0,x1) min(1,exp(x1 - x0));

ind = 1;
AR = 0;
for ii = 2:N
%     A = acc(x_chain(ind,1),x_cand(ii,1));
    A = acc2(lnw_v(ind,1),lnw_v(ii,1));
    if (A > u(ii,1))
        x_chain(ii,1) = x_cand(ii,1);
        ind = ind + 1;
        AR = AR + 1; 
    else
        x_chain(ii,1) = x_chain(ind,1);    
    end
end
AR = AR/(N-1);
% AC = autocorr(x_chain);

%% Figures 
figure(1)
hold on
plot(x_cand)
plot(x_chain,'r')
hold off

figure(2)
[f,x] = ksdensity(x_chain);
ind = (x >= 0 & x <= 1);

hold on
plot(x(ind),f(ind))
plot(x(ind),p(x(ind)),'r')
hold off
diff = mean((f(ind)-p(x(ind))).^2)

x_chain_whole = x_chain;
x_chain = x_chain(BurnIn+1:N,1);
figure(3)
[f,x] = ksdensity(x_chain);
ind = (x >= 0 & x <= 1);

hold on
plot(x(ind),f(ind))
plot(x(ind),p(x(ind)),'r')
hold off
diff = mean((f(ind)-p(x(ind))).^2)

%% Perturbation of the current state with a uniform innovation
clear all;
close all;

N = 15000;
BurnIn = 5000;

Eps = 0.5;
innov = 2*Eps*rand(N,1)-Eps;
ll = @(xx) -0.5*(log(2*pi) + xx.^2);

u = rand(N,1);
x_chain = zeros(N,1);
x_chain(1,1) = innov(1,1);
old_ll = ll(x_chain(1,1));

AR = 0;
for ii = 2:N
    x_cand = x_chain(ii-1,1) + innov(ii,1);
    new_ll = ll(x_cand);
    A = min(1,exp(new_ll - old_ll));
    if (A > u(ii,1))
        x_chain(ii,1) = x_cand;
        AR = AR + 1; 
    else
        x_chain(ii,1) = x_chain(ii-1,1);    
    end
end
AR = AR/(N-1);


figure(1)
[f,x] = ksdensity(x_chain);
ind = (x >= -5 & x <= 5);

hold on
plot(x(ind),f(ind))
plot(x(ind),exp(ll(x(ind))),'r')
hold off
diff = mean((f(ind)-exp(ll(x(ind)))).^2)


%% Pseudo marginal MCMC
% noisy measurement of loglikelihood ~ estimation
clear all;
close all;

N = 10000;
BurnIn = 5000;
N = N + BurnIn;

Eps = 0.5;
innov = 2*Eps*rand(N,1)-Eps;
% correct
ll = @(xx) -0.5*(log(2*pi) + xx.^2);
% unbiased noise
ll_noise1 = @(xx) log(exprnd(1)) - 0.5*(log(2*pi) + xx.^2);
% biased noise but state independent
ll_noise2 = @(xx) log(exprnd(2)) - 0.5*(log(2*pi) + xx.^2);
% biased state dependent noise but noise expectation state independent
ll_noise3 = @(xx) log(gamrnd(0.1+10*xx.^2,1./(0.1+10*xx.^2))) - 0.5*(log(2*pi) + xx.^2);
% biased state dependent noise with noise expectation state dependent
ll_noise4 = @(xx) log(exprnd(0.1+10*xx.^2)) - 0.5*(log(2*pi) + xx.^2);


[x_chain0, AR0] = MH_unif_innov(N,ll,innov);
[x_chain1, AR1] = MH_unif_innov(N,ll_noise1,innov);
[x_chain2, AR2] = MH_unif_innov(N,ll_noise2,innov);
[x_chain3, AR3] = MH_unif_innov(N,ll_noise3,innov);
[x_chain4, AR4] = MH_unif_innov(N,ll_noise4,innov);

x_chain0 = x_chain0(BurnIn+1:N,1);
x_chain1 = x_chain1(BurnIn+1:N,1);
x_chain2 = x_chain2(BurnIn+1:N,1);
x_chain3 = x_chain3(BurnIn+1:N,1);
x_chain4 = x_chain4(BurnIn+1:N,1);

figure(1)
subplot(1,2,1)
    qqplot(x_chain0)
subplot(1,2,2)
    [f,x] = ksdensity(x_chain0);
    hold on
    plot(x,f)
    plot(x,exp(ll(x)),'r')
    hold off


figure(2)
subplot(1,2,1)
    qqplot(x_chain1)
subplot(1,2,2)
    [f,x] = ksdensity(x_chain1);
    hold on
    plot(x,f)
    plot(x,exp(ll(x)),'r')
    hold off

figure(3)
subplot(1,2,1)
    qqplot(x_chain2)
subplot(1,2,2)
    [f,x] = ksdensity(x_chain2);
    hold on
    plot(x,f)
    plot(x,exp(ll(x)),'r')
    hold off
    
figure(4)
subplot(1,2,1)
    qqplot(x_chain3)
subplot(1,2,2)
    [f,x] = ksdensity(x_chain3);
    hold on
    plot(x,f)
    plot(x,exp(ll(x)),'r')
    hold off
    
figure(5)
subplot(1,2,1)
    qqplot(x_chain4)
subplot(1,2,2)
    [f,x] = ksdensity(x_chain4);
    hold on
    plot(x,f)
    plot(x,exp(ll(x)),'r')
    hold off
    
figure(6)
hold on
plot(1:N,x_chain0,'k')
plot(1:N,x_chain1,'b')
plot(1:N,x_chain2,'g')
plot(1:N,x_chain3,'m')
plot(1:N,x_chain4,'r')
hold off