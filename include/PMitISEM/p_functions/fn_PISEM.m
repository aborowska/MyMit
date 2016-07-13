function mit = fn_PISEM(theta, mit, w, cont, GamMat, X)

    [N,d] = size(theta);
    H = size(mit.p,2);       

    P = false; % partial
    if (nargin == 6)
        r = size(X,2); % X is (N)x(r) - for each draw a vector of length r
        if (r > 1)
            P = true;
        end
    end
    
    % auxiliary matrices
    df_mat = repmat(mit.df,N,1);
    w_mat = repmat(w,1,H);
    z = zeros(N,H);         % z <- tilde z <- E(z|theta) expected student t indicator
    z_wg = zeros(N,H);      % z_wg <- tilde(z/w) <- E(z/w|theta) weighted student t indicator
    xi = zeros(N,H);        % xi <- wg_ln <- tilde(log(w)) <- E(log(w)|theta) expected student t indicator
    delta = zeros(N,H);     % delta <- tilde(1/w) <- E(1/w) expected (inverse) IG draw
    rho = zeros(N,H);       % (theta-mu)'*Sigma^(-1)*(theta-mu)     
    
%% 1 step: EXPECTATION
% fprintf('Expectation \n')
    % update conditional expectation given last parameters
    for h = 1:H
        Sigma_h = mit.Sigma(h,:);
        df_h = mit.df(h);  
        mu_h = mit.mu(h,:);
        mit_h = struct('p',1,'mu',mu_h,'Sigma',Sigma_h,'df',df_h);
        if ~P
            z(:,h) = exp(log(mit.p(h)) + dmvgt(theta,mit_h,true, GamMat));
        else % if partial mitisiem then mit.mu is de fecto beta and mu_h is depends on X and is different for different draws 
%             mu_h = reshape(mit.mu(h,:),r,d); % this is beta
%             mu_h = X*mu_h; % (N)x(r) * (r)x(d) = (N)x(d) % this is a matrix of means
%             for ii = 1:N % MEX IT LATER TO AVIOD THE LOOP!!!
%                 mit_h = struct('p',1,'mu',mu_h(ii,:),'Sigma',Sigma_h,'df',df_h);
%                 z(ii,h) = exp(log(mit.p(h)) + dmvgt(theta(ii,:), mit_h, true, GamMat));
%             end
            z(:,h) = exp(log(mit.p(h)) + p_dmvgt(theta,mit_h,true, GamMat, X));
        end     
    end
    z = z./repmat(sum(z,2),1,H); % normalize indicator probability
    % psi(x) computes the digamma function of x
    % calculate digamma fucntions for eq. 11
    psi1 = psi((d+mit.df)/2);
    psi2 = psi(mit.df/2);
    for h = 1:H
        % calculate rho for eq. 10
        if ~P
            mu_h = mit.mu(h,:);
            mu_mat = repmat(mu_h,N,1);       % this is a matrix of means
        else
            mu_h = reshape(mit.mu(h,:),r,d); % this is beta
            mu_mat = X*mu_h;                 % (N)x(r) * (r)x(d) = (N)x(d) % this is a matrix of means    
        end
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
% fprintf('Maximisation\n')

    Sigma = zeros(H,d^2);
    if ~P
        mu = zeros(H,d);
    else
        mu = zeros(H,d*r);% this is beta (r)x(d), each row corresponds to each component and is reshaped to a (1)x(r*d) vector
    end
    
    for h = 1:H
        % update mu
        tmp_wg = w .* z_wg(:,h);
        tmp = repmat(tmp_wg,1,d);
        if ~P % standard mitisem
            mu(h,:) = sum(tmp.*theta,1) / sum(tmp_wg); 
            tmp_theta = theta - repmat(mu(h,:),N,1);
        else % partial mitisem

% CAN fn_beta BE USED??

            tmp_r = repmat(tmp_wg,1,r).*X; % (N)x(r)
            beta = X'*tmp_r; % (r)x(N) * (N)x(r) = (r)*(r) % <-- this is the 'denominator'
            if r == 2
                beta = [beta(2,2), -beta(1,2); -beta(2,1), beta(1,1)]/(beta(1,1)*beta(2,2)-beta(1,2)*beta(2,1));
                beta = beta*tmp_r'*theta; % here beta is already inverted
            else
                beta = beta\tmp_r'*theta; % inv(beta)=beta\ % (r)x(r) * (r)x(N) * (N)x(d) = (r)x(d) 
            end
            mu(h,:) = reshape(beta,1,r*d); % back to the mitisem format of storing parameters
            
            tmp_mu = X*beta; % mu = X*beta % <-- (N)x(r) * (r)x(d) = (N)x(d)
            tmp_theta = theta - tmp_mu;
        end
        tmp_Sigma = tmp_theta'*(tmp.*tmp_theta); % the denominator
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
            f_r = @(nu) -psi(nu/2) + log(nu/2) + 1 - w_exp(h);  
            if (f_r(tmp(2)) > 0) 
                df(h) = tmp(2);
            elseif  (f_r(tmp(1)) < 0) 
                df(h) = tmp(1);
            else
                df(h) = fzero(f_r, tmp);
            end
        end
    end
%%
%     if ~P % standard mitisem
        mit.mu = mu;
%     else % partial mitisem
%         mit.mu = reshape(beta,1,r*d);
%     end
    mit.Sigma = Sigma;
    mit.df = df;
    mit.p = eta;
end