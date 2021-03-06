function [theta_sim, lng_y, lnw, eps_bar, eps_sim, C_T] = NAIS_loglik_hl(y, par_SV, par_NAIS, cont, x_in)
%     n = size(y,1);
    b = par_NAIS.b;
    C = par_NAIS.C;
    [n, S] = size(b);
%     y = repmat(y,S);
    y = kron(y,ones(1,S));
    
    a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
%     a = 0.5*(log(abs(C))  - (b.^2)./C); 
    y_star = b./C;

    % Given the optimal IS parameters obtain the smoothed mean of signal for the importance model    
    par_KFS = IS_Model(par_NAIS, par_SV);
%     [~, ~, v1, F_inv1, eps_smooth1, K1, L1] = KFS(y_star,par_KFS); 
    [~, ~, v, F_inv, ~, K, L] = KFS_mex(y_star, par_KFS.P1, par_KFS.c, ...
                                         par_KFS.H, par_KFS.Q, par_KFS.d, par_KFS.T, par_KFS.R, par_KFS.Z);
     
    % LogLikelihood evaluation
    % (of the approximation linear state space model via KF)
    lng_y = -0.5*(n*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
%     lng_y = -0.5*(sum((v.^2).*F_inv) + sum(log(abs(F_inv))));
    
    % simulation smoothing for linear state space model to sample a signal trajectory via JSDK
	if (nargin == 4)
        m = size(par_SV,2); 
        % sample a signal trajectory via the simulation smoother
        if (((cont.err == 'n') && (m == 3)) || ((cont.err == 't') && (m == 4)))
            RND = randn(n,S);
            theta_sim = SimSmooth(y_star, v, F_inv, K, L, par_KFS, RND);
%         theta_sim = IS_sim(y_star, v, F_inv, K, L, par_KFS);
        else
            RND = randn(n,S);
            x_T = par_SV(:,end)';
            [theta_sim, eps_bar, eps_sim, C_T]= SimSmooth_hl(y_star, v, F_inv, K, L, par_KFS, RND);
            eps_bar = eps_bar';
            eps_sim =  eps_sim';
            C_T = C_T';
        end
    else
        theta_sim = x_in';
    end
    % compute the logweights
    if (cont.err == 'n')
    	lnp =  -0.5*(log(2*pi) +  theta_sim  + (y.^2)./exp(theta_sim));
%        	lnp =  -0.5*(theta_sim  + (y.^2)./exp(theta_sim));
    else % if (cont.err == 't')
    	nu = par_SV(:,4);
        nu = kron(nu',ones(n,1));
        p_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2); %% 
        y2 = (y.^2)./((nu-2).*exp(theta_sim));
        lnp = p_const - 0.5*(theta_sim + (nu+1).*log(1 + y2)); 
%         lnp = - 0.5*(theta_sim + (nu+1)*log(1 + y2));   
    end
    lng = a + b.*theta_sim - 0.5*C.*theta_sim.^2;
    lnp = sum(lnp,1);
    lng = sum(lng,1);
%     lnw = (lnp - lng)/(10*n);
    lnw = (lnp - lng);
    
%     lng_y = lng_y'/(10*n);
    lng_y = lng_y';
    lnw = lnw';
    if (nargin == 4)
        theta_sim = theta_sim';
    else
        theta_sim = x_in;
    end
    theta_sim = theta_sim(:,end); % pass only the last values for memory saving
end