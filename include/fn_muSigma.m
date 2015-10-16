function [mu, Sigma] = fn_muSigma(theta, w, mu)
% update mode and scale of the multivariate Student t matrix from IS weight
% theta - (T* x d) matrix of draws with highest IS weights
% w - vector of size T* of IS weights correponding to theta
% mu, Sigma - optimized mode and scale of the student t components    
    [N,d] = size(theta);
    w = exp(log(w) - log(sum(w))); % normalize weights
    tmp_w = repmat(w,1,d);
    if (nargin == 2) % if no mu supplied, compute mu
        mu = sum(tmp_w.*theta,1);
    end
    tmp_theta = theta - repmat(mu, N,1);
    
%     Sigma = upSigma(tmp_theta,w);
    Sigma = tmp_theta'*(tmp_w.*tmp_theta);

    % return scale in MitISEM format
    Sigma = reshape(Sigma,1,d^2);
end
