function d = posterior_debug(x, y, a, b)
end

function R = prior_debug(x, y, a, b)
% conjugate prior for sigma2 for normal likelihood with known mu = 0
    r1 = (x > 0);
    r2 = -Inf*ones(length(x),1);
    
    r2(r1==true) = a*log(b) - log(gamma(a)) - (a+1).*log(x) - b./x;
end