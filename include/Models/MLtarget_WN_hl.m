function [d, prior] = MLtarget_WN_hl(eps, sigma2_used, VaR_prelim)
    H = size(eps,2);
    
    prior =  prior_WN_hl(eps, sigma2_used, VaR_prelim); 
    d = prior(:,2) - 0.5*(H*log(2*pi) + H*log(sigma2_used) + sum((eps.^2)./sigma2_used,2));
end

function R = prior_WN_hl(eps, sigma2_used, VaR_prelim)
    N = size(eps,1);

    r1 = (fn_PL(sqrt(sigma2_used).*eps) <= VaR_prelim);   
    r2 = -Inf*ones(N,1);
    r2(r1==1,1) = 0;
    
    R = [r1, r2];
end