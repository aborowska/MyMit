function d = posterior_debug(x, y, a, b)
% log posterior for sigma2 for normal likelihood with known mu = 0
    N = size(x,1);
    T = size(y,1);
    prior =  prior_debug(x, a, b); 
    pdf = -Inf*ones(N,1);
    ind = (prior(:,1)==true);

%     pdf(ind) = - 0.5*T*log(x(ind,:)) -0.5*sum(y.^2)./x(ind,:);
    pdf(ind) = -0.5*(T*log(2*pi) + T*log(x(ind,:)) + sum(y.^2)./x(ind,:));

    d = prior(:,2) + pdf;
end

function R = prior_debug(x, a, b)
% conjugate prior for sigma2 for normal likelihood with known mu = 0
    r1 = (x > 0);
    r2 = -Inf*ones(length(x),1);
    
    ind = (r1 == true); 
    if (a == 0) % flat prior
        r2(ind,:) = 1; 
    else % conjugate prior inv gamma
        r2(ind,:) = a*log(b) - log(gamma(a)) - (a+1).*log(x(ind,:)) - b./x(ind,:);
    end
    R = [r1, r2];
end