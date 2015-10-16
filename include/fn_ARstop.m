function [AR_new, hstop] = fn_ARstop(w, AR_old, AR_tol)
% compute approximate expected acceptance rate from IS weights
% and indicator to finalize the number of mixture components in MitISEM
    AR_new = mean(w)/max(w);
    hstop = (abs((AR_new - AR_old)/AR_old) <= AR_tol);
end