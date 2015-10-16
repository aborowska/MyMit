function w = fn_ISwgts(lnk, lnd, norm)
% calculate IS weights from log (kernel) and (log) candidate density
% evaluations at N parameter (?) values 
% lnk - N vector of log kernel density evaluations
% lnd - N vector of log candidate density evaluations
% w - N vector of IS weights 
    w = lnk - lnd;
    w = w - max(w); % robustification
    w = exp(w);
    w(isnan(w)) = 0;
    if norm
        w = w/sum(w); % weights sum up to 1
    end
end