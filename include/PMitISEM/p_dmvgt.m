function dens = p_dmvgt(theta, mit, L, GamMat, X)
% density of a mixture of multivariate t distributions
% where the component means are varying between draws theta
% each mean is given by beta*X 
% where beta = reshape(mit.mu,r,d) 
% with: r = size(X,2)
%       d = size(theta,2) 
    L = double(L);
    dens = p_dmvgt_mex(theta, mit.mu, mit.Sigma, mit.df, mit.p, X, GamMat, L);    
end