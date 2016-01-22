function beta = fn_beta(theta, w, X)   
    r = size(X,2);
    w = exp(log(w) - log(sum(w))); % normalize weights

    tmp_r = repmat(w,1,r).*X; % (N)x(r)
    beta = X'*tmp_r; % (r)x(N) * (N)x(r) = (r)*(r) % <-- this is the 'denominator'
    beta = beta\tmp_r'*theta;    
end