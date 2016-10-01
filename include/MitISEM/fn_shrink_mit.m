function [crash, mit_s] = fn_shrink_mit(mit, tol_pr)
% shrinks a mixture of multivariate t densities (mit)     
% components are removed if their scale matrices are not pds (positive definite symmetric)
    H = length(mit.p);
    % indicator for zero variance for each component
    crash_Sigma = zeros(H,1); 
	for h = 1:H
        crash_Sigma(h,1) = fn_testSigma(mit.Sigma(h,:)); % indicator for crashing variance
    end
    % indicator for two small mixture probability
    crash_Prob = fn_testProb(mit.p, tol_pr);
    crash = (crash_Sigma | crash_Prob);

    if (all(crash))
        error('Mixture density shrinks to NULL, try different probability tolerance or increse N.');
    end
    
    if (any(crash)) % http://matlabtricks.com/post-23/tutorial-on-matrix-indexing-in-matlab
        % shrank mixture density
        d = sqrt(size(mit.Sigma,2));
        ind = logical(~crash);
%         ind_m = logical(ones(1,d));
%         ind_S = logical(ones(1,d^2));
%         
        mit_s.mu = mit.mu(ind,:);
        mit_s.Sigma = mit.Sigma(ind,:);      
        mit_s.df = mit.df(ind);
        mit_s.p = mit.p(ind);
        mit_s.p = mit_s.p/sum(mit_s.p);       
    else 
        mit_s = mit;
    end
    crash = any(crash);
end


function r = fn_testProb(p,tol)
    r = (p <= tol);
    r = r';
end