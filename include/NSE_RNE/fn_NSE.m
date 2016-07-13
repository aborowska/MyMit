function [NSE, RNE] = fn_NSE(y_opt, w_opt, G, arg) 
% [NSE_prob_IS, RNE_VaR_IS] = fn_NSE(y_opt, w_opt, G, VaR_IS) ;
% G = @(a,b) (fn_PL(a)<=b);
%     g = (PL(y_opt)<=VaR_IS);
    if (nargin > 3)
        g = G(y_opt,arg);
    else
        g = G(y_opt);
    end
    [N_g, d_g] = size(g);
    w_norm = w_opt/sum(w_opt);
    tmp_w = repmat(w_norm,1,d_g);
    ghat = sum(tmp_w.*g,1);
    tmp = g - repmat(ghat,N_g,1);
    NSE = sqrt(sum((tmp_w.^2).*(tmp.^2),1)); % NSE of the probability estimator
%     naive = sum(tmp_w.*(tmp.^2),1)/N_g;
%     RNE = naive./(NSE.^2);
% Book ch.7
    var_g = sum(tmp_w.*(g.^2),1) - ghat.^2;  % estimated variance of g 
    % (results from the function g itself and the posterior density)
    RNE = var_g./(N_g*NSE.^2);
% Geweke 1989
%     var_g = sum(tmp_w.*tmp.^2,1)/sum(w_opt);  % estimated variance of g 
%     % (results from the function g itself and the posterior density)
%     RNE = var_g./(N_g*NSE.^2);
end