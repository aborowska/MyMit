function d = posterior_binomial(n, y, theta)
    G = size(theta,1);   

    A = y + 1;
    B = n -y + 1;
    
%     logpdf_beta = @(x) - log(beta(A, B)) + (A-1).*log(x) + (B-1).*log(1-x); 
    prior = prior_binomial(theta);
    ind = (prior(:,1)==true);
    I = sum(ind);
    logpdf_beta = repmat(- log(beta(A, B)),I,1) +...
        repmat((A-1),I,1).*log(theta(ind,:)) +...
        repmat((B-1),I,1).*log(1-theta(ind,:)); 
    
    d = -Inf*ones(G,1);
    d(ind,1) = sum(logpdf_beta,2);
end

function R = prior_binomial(theta)
    [G, K] = size(theta);   
    r1 = all((theta >= 0) & (theta <= 1),2);
    r2 = -Inf*ones(G,1);
    r2(r1==true) = 1;

    R = [r1, r2];
end