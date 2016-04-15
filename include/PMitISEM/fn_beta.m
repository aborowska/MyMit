function [beta, Sigma] = fn_beta(theta, w, X)   
    r = size(X,2);
    d = size(theta,2);
%     w = exp(log(w) - log(sum(w))); % normalize weights

    tmp_r = repmat(w,1,r).*X; % (N)x(r)
    beta = X'*tmp_r; % (r)x(N) * (N)x(r) = (r)*(r) % <-- this is the 'denominator'
    beta = beta\tmp_r'*theta;    
 
    tmp_mu = X*beta; % mu = X*beta % <-- (N)x(r) * (r)x(d) = (N)x(d)
    tmp_theta = theta - tmp_mu;

    tmp_r = repmat(w,1,d); % (N)x(r)
    tmp_Sigma = tmp_theta'*(tmp_r.*tmp_theta);
    tmp_Sigma = tmp_Sigma/sum(w);
    
    beta = reshape(beta,1,r*d);
    Sigma  = reshape(tmp_Sigma,1,d^2);
end