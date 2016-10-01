function d = posterior_binomial_hl(n, y, theta, PVR_prelim)
    G = size(theta,1);   

    A = y + 1;
    B = n -y + 1;
    
%     logpdf_beta = @(x) - log(beta(A, B)) + (A-1).*log(x) + (B-1).*log(1-x); 
    prior = prior_binomial_hl(theta, PVR_prelim);
    ind = (prior(:,1)==true);
    I = sum(ind);
    logpdf_beta = repmat(- log(beta(A, B)),I,1) +...
        repmat((A-1),I,1).*log(theta(ind,:)) +...
        repmat((B-1),I,1).*log(1-theta(ind,:)); 
    
    d = -Inf*ones(G,1);
    d(ind,1) = sum(logpdf_beta,2);
end

function R = prior_binomial_hl(theta, PVR_prelim)
    [G, ~] = size(theta);   
    c1 = all((theta >= 0) & (theta <= 1),2);
    c2 = (regret_binomial(theta) >= PVR_prelim);
    r1 = c1 & c2;
    r2 = -Inf*ones(G,1);
    r2(r1==true) = 1;

    R = [r1, r2];
end