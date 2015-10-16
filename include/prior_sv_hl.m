function r2 = prior_sv_hl(theta, prior_const)
    c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1);
    phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2);
    eta = theta(:,4);   % prior: eta ~ N(0,1)
    eps = theta(:,5);   % prior: eps ~ N(0,1)
    
    logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
    logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
    % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
    logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
    % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) - 0.5*log(x) - 0.5*x;    

%     c1 = (PL <= VaR);
    c2 = ((phi > 0) & (phi < 1));
%     c3 = (s2 > 0);
    c3 = (s2 > 5e-3);   
     
%     r1 = (c1 & c2 & c3);
    r1 = (c2 & c3);
    r2 = -Inf*ones(length(eta),1);
    
    r2(r1==true) = logpdf_norm(c(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
%     r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_invgamma(s2(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));
    r2(r1==true) = r2(r1==true) + prior_const(1,1) - 0.5*(eta(r1==true)).^2;
    r2(r1==true) = r2(r1==true) + prior_const(1,1) - 0.5*(eps(r1==true)).^2;

%     R = [r1, r2];
end