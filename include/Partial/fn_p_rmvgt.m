function [theta, lnk] = fn_p_rmvgt(N, pmit, d, partition, kernel, y, X)
% sampling from mixture of multivariate t densities
% with the parameter vector of size d divided into ordered subsets
% the mean may depend on the matrix X being a function of data y 
% and draws of parameters from previous subsets
 
    S = size(pmit,2);
    theta = zeros(N,d);
    
    for s = 1:S
        % the dimension of the s-th set
%         d_s = size(pmit(s).Sigma,2);
%         d_s = sqrt(d_s);
        [s1, s2] = fn_partition_ends(partition, d, s);
       
        if (s == 1)
            % Mixtue of t (mit) parameters:
            mu = pmit(s).mu;
            Sigma = pmit(s).Sigma;
            df = pmit(s).df;
            p = pmit(s).p;
            theta(:,s1:s2) = rmvgt2(N, mu, Sigma, df, p); % Sampling from the mixture of t
        else
            X = fn_const_X(theta(:,1:s1-1));  % X is a function of the draws from previous subsets
            beta = pmit(s).mu;
            Sigma = pmit(s).Sigma;
            df = pmit(s).df;
            p = pmit(s).p;
            theta(:,s1:s2) = p_rmvgt2(N, beta, X, Sigma, df, p); % Sampling from the mixture of t
        end
    end
    
    % lnk - N vector of log-kernel evaluations at draws
    fprintf('\n'); 
    fprintf('kernel computation')
    fprintf('\n'); 
    lnk = kernel(theta); 
end 