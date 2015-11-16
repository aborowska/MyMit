function d = tail_normal(x,VaR_prelim)
    N = size(x,1);
    d = -Inf*ones(N,1);
    r1 = (x<VaR_prelim);
    d(r1==1) = - 0.5*(log(2*pi) + x(r1==1).^2);
end