function [theta, lnd, input_X] = fn_p_rmvgt_dpmit3(N, pmit, d, S, partition, fn_const_X, fn_input_X, GamMat)
% sampling from mixture of multivariate t densities
% with the parameter vector of size d divided into ordered subsets
% the mean may depend on the matrix X being a function of data y 
% and draws of parameters from previous subsets
% AND
% density evaluation of theta on the Partial MitISEM structure
% theta is partitioned according to partition
    theta = zeros(N,d);

    % Sample from the partial candidate and evaluate on it
    for s = 1:S
        [s1, s2] = fn_partition_ends(partition, d, s);
       
        if (s == 1)
            % Mixtue of t (mit) parameters:
            mu = pmit(s).mu;
            Sigma = pmit(s).Sigma;
            df = pmit(s).df;
            p = pmit(s).p;
            theta(:,s1:s2) = rmvgt2(N, mu, Sigma, df, p); % Sampling from the mixture of t
            input_X = fn_input_X(theta(:,s1:s2));
            % evaluate on the candidate
            lnd = dmvgt(theta(:,s1:s2), pmit(s), true, GamMat);
        else
            input_X.theta = theta(:,1:s1-1);       
            [X, input_X] = fn_const_X(input_X);  % X is a function of the draws from previous subsets
            beta = pmit(s).mu;
            Sigma = pmit(s).Sigma;
            df = pmit(s).df;
            p = pmit(s).p;
            theta(:,s1:s2) = p_rmvgt2(N, beta, X, Sigma, df, p); % Sampling from the mixture of t
            % evaluate on the candidate
            lnd = lnd + p_dmvgt(theta(:,s1:s2), pmit(s), true, GamMat, X);
        end
    end
    input_X.theta = theta;
    [~, input_X] = fn_const_X(input_X);
end