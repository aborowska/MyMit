function d = MLtarget_arch_hl(alpha, eps, y_T, S, VaR, L)
    [N,H] = size(eps);
    alpha = repmat(alpha,N,1);
    
    if (H == 1)
        h = S*(1-alpha) + (y_T^2)*alpha;
        c1 = (fn_PL(sqrt(h).*eps) <= VaR);
    else
        y_T1 = predict_arch(alpha, y_T, S, H, eps);
        c1 = (fn_PL(y_T1) <= VaR);
    end
    
    d = -Inf*ones(N,1);
    for ii = 1:N
        if c1(ii,1)
            d(ii,1) = - sum(0.5*(log(2*pi) + eps(ii,:).^2),2);
        end
    end
    
    if (~L)
        d = exp(d - max(d));
    end
end
 