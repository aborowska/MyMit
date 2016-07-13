function d = posterior_debug(sigma2, y, a, b, L)
% log posterior for sigma2 for normal likelihood with known mu = 0
    sigma2 = sigma2(:,1);
    N = size(sigma2,1);
    T = size(y,1);
    prior =  prior_debug(sigma2, a, b); 
    pdf = -Inf*ones(N,1);
    ind = (prior(:,1)==true);

%     pdf(ind) = - 0.5*T*log(x(ind,:)) -0.5*sum(y.^2)./x(ind,:);
    pdf(ind) = -0.5*(T*log(2*pi) + T*log(sigma2(ind,:)) + sum(y.^2)./sigma2(ind,:));

    d = prior(:,2) + pdf;
    if (~L)
        d = exp(d-max(d));
    end
end

function R = prior_debug(sigma2, a, b)
% conjugate prior for sigma2 for normal likelihood with known mu = 0
    r1 = (sigma2 > 0);
    r2 = -Inf*ones(length(sigma2),1);
    
    ind = (r1 == true); 
    if (a == 0) % flat prior
        r2(ind,:) = 1; 
    else % conjugate prior inv gamma
        r2(ind,:) = a*log(b) - log(gamma(a)) - (a+1).*log(sigma2(ind,:)) - b./sigma2(ind,:);
    end
    R = [r1, r2];
end