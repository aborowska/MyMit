function dens = p_dmvgt(theta, mit, L, GamMat, X)

    L = double(L);
    dens = p_dmvgt_mex(theta, mit.mu, mit.Sigma, mit.df, mit.p, X, GamMat, L);    
end