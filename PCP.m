N = 10000;
sigma1 = 1;
sigma2 = 2;
% sample membership
memb = randsample(0:1,N,true,[0.5 0.5])';
eps = abs(randn(N,1));
eps1 = sigma1.*eps;
eps2 = -sigma2.*eps;
hist([eps1,eps2])
eps = memb.*eps1 + (1-memb).*eps2;
hist(eps,100)

% simple AR(1)
rho = 0.8;
y = zeros(N,1);
y(1,1) = eps(1,1);
for ii = 2:N
    y(ii,1) = rho*y(ii-1,1) + eps(ii);
end