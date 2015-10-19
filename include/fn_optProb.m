function p = fn_optProb(p, lnk, lnd, cont)
    p_nc = cont.pnc;
%     lambda0 = log((1-p_nc).*p./p_nc);
    p = [(1-p_nc).*p, p_nc];                % p=[0.9 0.1] % p = [(1-p_nc).*mit.p, p_nc];
    
%     lambda0 = log(p(1:(H-1))/p(H));
   % log(fn_lambda(lambda0))
%     fn_lambda = @(x) [exp(x),1]./sum(exp(x)+1);
  
%     fn_lnf = @(x)  fn_CoVsq(log(fn_lambda(x)), lnk, lnd);

%   fn_lnf = @(x)  fn_CoVsq_C(log(x),lnk_up(1:Np,:), lnd_up);
    fn_lnf = @(x)  fn_CoVsq_C(log(x), lnk, lnd);
    
%     function [f, g] = fn_lnf(x) %, lnk, lnd)
%         [f, g] = fn_CoVsq(log(fn_lambda(x)), lnk, lnd);
%     end

    options = optimset('GradObj','on','Display','off');
    try 
        p = fminunc(fn_lnf, p, options);
%     p = fminunc(fn_lnf, lambda0, options);
    end
%     p = fn_lambda(p);
    p = exp(p)/sum(exp(p));
end

% tic
% p_nc = cont.pnc;
% p = [(1-p_nc).*p, p_nc];                % p=[0.9 0.1]
% lambda0 = log(p) - log(1-p);
% fn_lnf = @(x)  fn_CoVsq_C(x, lnk, lnd); % fn_lnf = @(x) fn_CoVsq_C(x, lnk_up(1:Np,:), lnd_up);
% options = optimset('GradObj','on');
% p = fminunc(fn_lnf, lambda0, options);
% p = exp(p)./(1+exp(p));
% p = p./sum(p);
% time = toc;
% end