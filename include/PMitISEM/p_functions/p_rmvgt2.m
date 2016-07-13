function theta = p_rmvgt2(N, beta, X, Sigma, df, p)
% Random sampling from mixture of t densitites with diffente means 
% each draw has mean mu^i = beta*X^i
% N - number of draws
% mit - list with parameters of mixture of t density
    [H, d] = size(Sigma); % number of components, dimension of t distribution
    d = sqrt(d);
    r = size(X, 2); 
 
    % sample membership
    memb = randsample(1:H,N,true,p);
    % randsample(1:3,10,true,[0.1 0.3 0.6])
    theta = zeros(N, d);
    
    for h = 1:H
        ind_h = (memb == h);
        n_h = sum(ind_h);
        if (n_h>0)
            beta_h = reshape(beta(h,:),r,d);
            mu_h = X(ind_h,:)*beta_h;
            Sigma_h = Sigma(h,:);
            Sigma_h = reshape(Sigma_h,d,d);
            df_h = df(h); 
            draw_h = rmvt(mu_h,Sigma_h,df_h,n_h); 
            theta(ind_h,:) = draw_h;
        end
    end
end
