s2 = theta(:,3);
[~,ind] = sort(s2);

figure(1)
subplot(3,1,1)
plot(s2(ind))
subplot(3,1,2)
plot( logpdf_invgamma(s2(ind)))
subplot(3,1,3)
plot(s2(ind),exp(logpdf_invgamma(s2(ind))))

% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -0.5*log(2), -log(gamma(0.5))];

logpdf_gamma = @(a,b,x) -a*log(b) -log(gamma(a)) + (a-1)*log(x) - x/b;
logpdf_invgamma = @(a,b,x) a*log(b) - log(gamma(a)) - (a+1)*log(x) - b./x;
% logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;

a = 2.5;
d_gamma1 =  logpdf_gamma(a,0.025,1./s2);
d_gamma2 =  logpdf_gamma(a,40,1./s2);
d_invgam = logpdf_invgamma(a,0.025,s2);
% r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));

