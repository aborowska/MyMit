function [CV_new, hstop] = fn_CVstop(w, CV_old, CV_tol)
% compute new CV from IS weights
% and indicator to finalize the number of mixture components in MitISEM
    if std(w) == 0 
        error('IS weights w are constnat, try increasing number of draws N.')
    end
    CV_new = std(w)/mean(w);
    if ((nargin > 1) && (nargout > 1))
        hstop = (abs((CV_new - CV_old)/CV_old) <= CV_tol);
    end
end