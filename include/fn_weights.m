function w = fn_weights(theta, kernel, mit, L, GamMat)
    lnk = kernel(theta); % log-kernel evaluations at draws 
    lnd = dmvgt(theta, mit, true, GamMat); % density of mixture of t's at draws
    if L % log = true --> log weights
        w = lnk - lnd;
        w(isnan(w)) = -Inf;
    else % exponential weights
        w = fn_ISwgts(lnk, lnd, false); % false = not normalized 
    end
end