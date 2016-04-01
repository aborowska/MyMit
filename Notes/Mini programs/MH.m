clear all;
close all;
% p = @(xx) 6*(xx-xx.^2);
% q = @(xx) 1;
p = @(xx) exp((-xx.^2)/2)/sqrt(2*pi);
q = @(xx) 1;
f = @(xx) tan(pi*xx - pi/2);
N = 15000;
BurnIn = 5000;

x_cand = f(rand(N,1));

u = rand(N,1);
x_chain = zeros(N,1);
x_chain(1,1) = x_cand(1,1);

acc = @(x0,x1) min(1,exp(log(p(x1)) - log(p(x0)) + log(q(x1)) - log(q(x0))));
% p_val = p(x_cand);
% q_val = q(x_cand);

ind = 1;
AR = 0;
for ii = 2:N
    A = acc(x_chain(ind,1),x_cand(ii,1));
    if (A > u(ii,1))
        x_chain(ii,1) = x_cand(ii,1);
        ind = ind + 1;
        AR = AR+1; 
    else
        x_chain(ii,1) = x_chain(ind,1);    
    end
end
AR = AR/(N-1);
AC = autocorr(x_chain);

figure(1)
hold on
plot(x_cand)
plot(x_chain,'r')
hold off

figure(2)
[f,x] = ksdensity(x_chain);
% ind = (x >= 0 & x <= 1);
ind = (x >= -5 & x <= 5);

hold on
plot(x(ind),f(ind))
plot(x(ind),p(x(ind)),'r')
hold off
diff = mean((f(ind)-p(x(ind))).^2)

x_chain_whole = x_chain;
x_chain = x_chain(BurnIn+1:N,1);
figure(3)
[f,x] = ksdensity(x_chain);
% ind = (x >= 0 & x <= 1);
ind = (x >= -5 & x <= 5);

hold on
plot(x(ind),f(ind))
plot(x(ind),p(x(ind)),'r')
hold off
diff = mean((f(ind)-p(x(ind))).^2)
