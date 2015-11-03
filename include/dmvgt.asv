function dens = dmvgt(theta, mit, L, GamMat)
% density of a mixture of multivariate t distributions
% L (log) - return log-density values if L=true
    [H,d] = size(mit.mu); % number of components, dimension of t distribution
    [N,~] = size(theta);
    dcoms = zeros(N,H);
    for h = 1:H
        mu_h = mit.mu(h,:);
        Sigma_h = mit.Sigma(h,:);
        Sigma_h = reshape(Sigma_h,d,d);
        df_h = mit.df(h);
        for ii = 1:N
            dcoms(ii,h) =  dmvt(theta(ii,:), mu_h, Sigma_h, df_h, GamMat);
        end
    end
    tmp = log(repmat(mit.p,N,1)) + log(dcoms);
    dens = sum(exp(tmp),2);
    if (L == true)
        dens = log(dens);
    end
    
    L = double(L);
    lnd_mex = dmvgt_mex(theta, mit.mu, mit.Sigma, mit.df, mit.p, GamMat, 1);
    fprintf('\n *** sum(abs(lnd_mex-lnd)> eps)  = %6.4f ***\n\n', sum(abs(lnd_mex-dens)>eps) );

    
end