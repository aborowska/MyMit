function d = posterior_debug_hl(theta, y, a, b, VaR_prelim, L)
% log posterior for sigma2 for normal likelihood with known mu = 0
    [N, H] = size(theta);
    H = H-1; % number of epsilons i.e. the forecast horizon 
    T = size(y,1);
    sigma2 = theta(:,1);
    eps = theta(:,2:H+1);
    
    prior =  prior_debug_hl(sigma2, eps, a, b, VaR_prelim); 
    pdf = -Inf*ones(N,1);
    ind = (prior(:,1)==true);

%     pdf(ind) = - 0.5*T*log(x(ind,:)) -0.5*sum(y.^2)./x(ind,:);
    pdf(ind) = -0.5*(T*log(2*pi) + T*log(sigma2(ind,:)) + sum(y.^2)./sigma2(ind,:));

    d = prior(:,2) + pdf - 0.5*(H*log(2*pi) + sum(eps.^2,2));
    if (~L)
        d = exp(d-max(d));
    end
end

function R = prior_debug_hl(sigma2, eps, a, b, VaR_prelim)
% conjugate prior for sigma2 for normal likelihood with known mu = 0
    H = size(eps,2);
    c1 = (sigma2 > 0);
    c2 = (fn_PL(repmat(sqrt(sigma2),1,H).*eps) <= VaR_prelim);
    
    r1 = c1 & c2;
    r2 = -Inf*ones(length(sigma2),1);
    
    ind = (r1 == true); 
    if (a == 0)
        r2(ind,:) = 1;
    else
        r2(ind,:) = a*log(b) - log(gamma(a)) - (a+1).*log(sigma2(ind,:)) - b./sigma2(ind,:);
    end
    R = [r1, r2];
end