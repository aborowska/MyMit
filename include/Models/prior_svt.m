function r2 = prior_svt(theta, prior_const)
%     % prior is an Nx2 matrix: 
%     % 1 col - constraint satisfied?
%     % 2 col - prior val an the corresponding point
       
    c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1)
    phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5)
    s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2)
    nu = theta(:,4);    % prior: nu-2 ~ exp(1)
    
logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) - 0.5*log(x) - 0.5*x;    
    
    c1 = ((phi > 0) & (phi < 1));
    c2 = (s2 > 5e-3);
    c3 = (nu > 2);
    r1 = (c1 & c2 & c3);
    r2 = -Inf*ones(length(c),1);
    
    r2(r1==true) = logpdf_norm(c(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
%     r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_invgamma(s2(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));
    r2(r1==true) = r2(r1==true) - ( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  
    
%     R = [r1, r2];
end