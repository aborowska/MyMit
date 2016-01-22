function [CV_new, hstop] = fn_CVstop(w, CV_old, CV_tol, lnk, lnd, lnd_new)
% compute new CV from IS weights
% and indicator to finalize the number of mixture components in MitISEM
    if std(w) == 0 
        error('IS weights w are constnat, try increasing number of draws N.')
    end
    if (nargin > 3)
        w = lnk - lnd;
        M = max(w);       
        w = exp(w - M);
        w(lnk == -Inf) = 0;
        mean_i = mean(w);
%         mean_i = exp(M + log(mean_i));
        
        w = 2*lnk - lnd - lnd_new;
%         M = max(w);
        w = exp(w - 2*M);
        w(lnk == -Inf) = 0;
        var_i = mean(w);
%         var_i = exp(M + log(var_i));
        
        var_i = var_i - mean_i^2;
        CV_new = sqrt(var_i)/mean_i;
    else
        CV_new = std(w)/mean(w);
    end
    if ((nargin > 1) && (nargout > 1))
        hstop = (abs((CV_new - CV_old)/CV_old) <= CV_tol);
    end
end