function mit = fn_ISEM(theta, mit, w, cont, GamMat)

    [N,d] = size(theta);
    H = size(mit.p,2);
    
    % auxiliary matrices
    df_mat = repmat(mit.df,N,1);
    w_mat = repmat(w,1,H);
    z = zeros(N,H);         % z <- tilde z <- E(z|theta) expected student t indicator
    z_wg = zeros(N,H);      % z_wg <- tilde(z/w) <- E(z/w|theta) weighted student t indicator
    xi = zeros(N,H);        % xi <- wg_ln <- tilde(log(w)) <- E(log(w)|theta) expected student t indicator
    delta = zeros(N,H);     % delta <- tilde(1/w) <- E(1/w) expected (inverse) IG draw
    rho = zeros(N,H);       % (theta-mu)'*Sigma^(-1)*(theta-mu)     
    
%% 1 step: EXPECTATION
    % update conditional expectation given last parameters
    for h=1:H
        mu_h = mit.mu(h,:);
        Sigma_h = mit.Sigma(h,:);
        df_h = mit.df(h);
        mit_h = struct('p',1,'mu',mu_h,'Sigma',Sigma_h,'df',df_h);
        z(:,h) = exp(log(mit.p(h)) + dmvgt(theta,mit_h,true, GamMat));
    end
    z = z./repmat(sum(z,2),1,H); % normalize indicator probability
    % psi(x) computes the digamma function of x
    % calculate digamma fucntions for eq. 11
    psi1 = psi((d+mit.df)/2);
    psi2 = psi(mit.df/2);
    for h=1:H
        % calculate rho for eq. 10
        mu_h = mit.mu(h,:);
        mu_mat = repmat(mu_h,N,1);
        Sigma_h = mit.Sigma(h,:);
        Sigma_h = reshape(Sigma_h,d,d);
        tmp = chol(inv(Sigma_h));
        tmp = tmp*(theta - mu_mat)';
        rho(:,h) = sum(tmp.^2,1);
        % calculate xi  E(log(w)|theta) for eq. 11
        tmp_z = [z(:,h),1-z(:,h)];
        tmp_l = [log((rho(:,h) + df_mat(:,h))/2), log(df_mat(:,h)/2)];
        tmp_psi = repmat([psi1(h),psi2(h)],N,1);
        tmp_l = tmp_l - tmp_psi;
        xi(:,h) = sum(tmp_l .* tmp_z,2);  
    end
    % calculate weighted membership eq. 10
    z_wg = z .* (d + df_mat)./(rho + df_mat);
    % calculate delta eq. 12 using equality in eq. 10 
    delta = z_wg + (1-z);
    
%% 2 step: MAXIMIZATION
    Sigma = zeros(H,d^2);
    mu = zeros(H,d);
    for h=1:H
        % update mu
        tmp_wg = w .* z_wg(:,h);
        tmp = repmat(tmp_wg,1,d);
        mu(h,:) = sum(tmp.*theta,1) / sum(w .* z_wg(:,h));
        tmp_theta = theta - repmat(mu(h,:),N,1);
        % update sigma
%         tmp_Sigma = arrayfun(@(ii) tmp_theta(ii,:)'*tmp_theta(ii,:), 1:N, 'un', 0);
%         tmp_Sigma = cat(3,tmp_Sigma{:});
%         tmp = zeros(1,1,N); tmp(1,1,:) = tmp_wg; tmp = repmat(tmp,d,d,1); 
%         tmp_Sigma = sum(tmp .* tmp_Sigma,3)/ sum(w .* z(:,h));
%         Sigma(h,:) = reshape(tmp_Sigma,1,4);


%          tmp_Sigma = upSigma(tmp_theta,tmp_wg)/ sum(w .* z(:,h));
        tmp_Sigma = tmp_theta'*(tmp.*tmp_theta);
        tmp_Sigma = tmp_Sigma/ sum(w .* z(:,h));
    
         
         
         
         Sigma(h,:) = reshape(tmp_Sigma,1,d^2);
    end
    % eta = updated probability (of theta^i belonging to the h-th component)   
    eta = sum(w_mat.* z,1)/ sum(w); 
    
%% update degrees of freedom
    df = mit.df; % if not optimized df are not alterned
    if cont.df.opt
        % calculate weighted expectation in eq. 16
        w_exp = sum(w_mat.*(xi + delta),1)/sum(w);
        tmp = cont.df.range;
        for h=1:H  
            % decreasing objective function
            r = @(nu) -psi(nu/2) + log(nu/2) + 1 - w_exp(h);  
            if (r(tmp(2)) > 0) 
                df(h) = tmp(2);
            elseif  (r(tmp(1)) < 0) 
                df(h) = tmp(1);
            else
                df(h) = fzero(r, tmp);
            end
        end
    end
%%
    mit.mu = mu;
    mit.Sigma = Sigma;
    mit.df = df;
    mit.p = eta;
end