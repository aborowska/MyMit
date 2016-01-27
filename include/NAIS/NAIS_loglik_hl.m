function [theta_sim, lng_y, lnw, eps_bar, eps_sim, C_T, lnp_T, RND] = NAIS_loglik_hl(y, par_SV, par_NAIS, cont)
% Computes the marginal likelihood from the KF and the NAIS importance
% weight for the signal generated from the simulation smoother.
% If x_T are passed in - substitute the last draws of state from the simulation
% smoother with them.
% For the draw of the last state, its log evalaution on the Gaussian
% candidate is passed as output.
    m = size(par_SV,2);
    b = par_NAIS.b;
    C = par_NAIS.C;
    [n, S] = size(b);
    y = kron(y,ones(1,S));

    a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
    y_star = b./C;

    % Given the optimal IS parameters obtain the smoothed mean of signal for the importance model    
    par_KFS = IS_Model(par_NAIS, par_SV);
    [~, ~, v, F_inv, ~, K, L] = KFS_mex(y_star, par_KFS.P1, par_KFS.c, ...
                                         par_KFS.H, par_KFS.Q, par_KFS.d, par_KFS.T, par_KFS.R, par_KFS.Z);
     
    % LogLikelihood evaluation
    % (of the approximation linear state space model via KF)
    lng_y = -0.5*(n*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
    
    % simulation smoothing for linear state space model to sample a signal trajectory via JSDK
    RND = randn(n,S);
%     if (((cont.err == 'n') && (m==6)) || (((cont.err == 't')) && (m==7)))
%         RND(n,:) = par_SV(:,end)'; 
%     end
    [theta_sim, eps_bar, eps_sim, C_T]= SimSmooth(y_star, v, F_inv, K, L, par_KFS, RND);

    
    if (((cont.err == 'n') && (m==6)) || (((cont.err == 't')) && (m==7)))
        theta_sim(n,:) = par_SV(:,end)'; 
    end
    
% transition probability    
mu_T = theta_sim(n,:) - par_SV(:,1)'-(par_SV(:,2)').*(theta_sim(n-1,:)- par_SV(:,1)');
lnp_T = -0.5*(log(2*pi) + log(par_SV(:,3)') + (mu_T.^2)./(par_SV(:,3)')); 
lnp_T = lnp_T';

    % compute the logweights
    if (cont.err == 'n')
    	lnp =  -0.5*(log(2*pi) +  theta_sim  + (y.^2)./exp(theta_sim));
    else % if (cont.err == 't')
    	nu = par_SV(:,4);
        nu = kron(nu',ones(n,1));
        p_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2);  
        y2 = (y.^2)./((nu-2).*exp(theta_sim));
        lnp = p_const - 0.5*(theta_sim + (nu+1).*log(1 + y2)); 
    end
    lng = a + b.*theta_sim - 0.5*C.*theta_sim.^2;
%      lnp_T = lnp(n,:)';
%     lng_T = lng(n,:)';
    lnp = sum(lnp,1);
    lng = sum(lng,1);
    lnw = (lnp - lng);
    
    % transpose for compliance with storing of draws in QERMit
    lng_y = lng_y';
    lnw = lnw'; 
    eps_bar = eps_bar';
    eps_sim =  eps_sim';
    C_T = C_T';
    RND = RND(end,:)';
    theta_sim = theta_sim((end-cont.HP):end,:)';  % pass only the last values for memory saving
end