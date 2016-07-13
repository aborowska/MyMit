function [ghat, NSE, RNE, naive] = Mit_IS(N, kernel, mit, G)
% inpout: 
% perforem IS using mixture of t as the imnpostance density
% kernel - function which computes the kernel
% mit - mixture parameters
% G - function used in importnace sampling, e.g. identity: G(theta)=theta
% --> mean IS estiamtor
% output:
% ghat - IS estimator
% NSE - numerical standard error estimator
% RNE - relative numerical efficiency estiamtor
    
    [theta, lnk, ~, ~] = fn_rmvgt_robust(N, mit, kernel);
    lnd = dmvgt(theta,mit,true);
    norm = false;
    w = fn_ISwgts(lnk, lnd, norm);
%     [N, d] = size(theta);
    
    g = G(theta);   
	[N, d] = size(g);
    
    w_norm = w/sum(w);
    tmp_w = repmat(w_norm,1,d);
    ghat = sum(tmp_w.*g,1);
    tmp = g - repmat(ghat,N,1);
    NSE = sqrt(sum((tmp_w.^2).*(tmp.^2),1));
    naive = sum(tmp_w.*(tmp.^2),1)/N;
    RNE = naive./(NSE.^2);
    %%%%
%     tmp_w = repmat(w,1,d);
%     ghat = sum(tmp_w.*g,1)/sum(w);
%     tmp = g - repmat(ghat,N,1);
%     NSE = sqrt(sum((tmp_w.*tmp).^2,1));
%     naive = (sum(tmp_w.*tmp.^2,1)/N );
%     RNE = naive ./ (NSE.^2);  % naive variance/NSE^2
    % if target and importance density coincide, RNE equals one, if the
    % candidate is very poor, RNE is close to zero
end