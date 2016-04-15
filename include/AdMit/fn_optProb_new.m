function p_opt = fn_optProb_new(p, lnk, lnd, cont)
    p_nc = cont.pnc;
    p = [(1-p_nc).*p, p_nc];                % p=[0.9 0.1] % p = [(1-p_nc).*mit.p, p_nc];
    
    % unconstained optimisation --> optimise wrt the transformed parameters
    % transformation to have positivity: exp
    % transformation to have summability: normalisation
    fn_trans  = @(x) exp(x)/sum(exp(x)); 
    fn_trans_back = @(x) log(x);
%     p_trans = log(p/p(end));  % normalise wrt to the last parameter
    p_trans = fn_trans_back(p);
    
    
    fn_trans(fn_trans_back(p));
    fn_lnf = @(x) fn_CoVsq_C_new(x, lnk, lnd); % log because we work withlog-probabilies
% fn_lnf = @(x) fn_CoVsq_C_new(x, lnk_up(1:Np,:), lnd_up);
  
    options = optimset('GradObj','on','Display','off');    
    p_opt_trans = fminunc(fn_lnf, p_trans, options);
    p_opt = fn_trans(p_opt_trans);
  
end