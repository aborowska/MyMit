function [CV_new, AR_new, hstop] = fn_Mitstop(meth, w , CV_old, CV_tol, AR_old, AR_tol)
% indicates end of adding mixture components in MitISEM candidate
% assesses convergence of the MitISEM algorithm depending on the method (CV/AR)
% meth - the chosen convergence criterion, CV or AR
    meth_set = {'CV','AR'};
    if ~ismember(meth, meth_set)
        error('Stopping method for MitISEM should be CV or AR.')
    end
    
    [CV_new, CV_stop] = fn_CVstop(w, CV_old, CV_tol);
	[AR_new, AR_stop] = fn_ARstop(w, AR_old, AR_tol);
    
    if  strcmp(meth,'CV') %  compares strings
        hstop = CV_stop;
    else
        hstop = AR_stop;
    end
end

