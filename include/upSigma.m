function Sigma = upSigma(theta, w)
% Weighted covariance estimate, i.e. sum_i=0^N (w * theta_i^T * theta_i) 
    [N,d] = size(theta);
    Sigma = arrayfun(@(ii) theta(ii,:)'*theta(ii,:), 1:N, 'un', 0);
    Sigma = cat(3,Sigma{:});
    
    tmp = zeros(1,1,N); 
    tmp(1,1,:) = w; 
    tmp = repmat(tmp,[d,d,1]); 
    
    Sigma = sum(tmp .* Sigma,3);
end
