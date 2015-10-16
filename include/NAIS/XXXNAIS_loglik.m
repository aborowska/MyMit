function [lnL_hat, x, lnw_x, a, b, C] = NAIS_loglik(par_SV, par_NAIS, y, cont, RND)
    % Algorithm 1: IS using an approximation linear state model
    % (with antithetic variables for variance reduction)
    n = length(y);
    S = cont.S;
    fprintf('%6.4f\t',par_SV); 
    
    % Algorithm 2: Efficient importance parameters via NAIS
    par_NAIS = NAIS_param(par_NAIS, y, par_SV, cont); 

    b = par_NAIS.b;
    C = par_NAIS.C;
    a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
%     a = 0.5*(log(abs(C))  - (b.^2)./C); <-- was

    y_star = b./C;

    % Given the optimal IS parameters obtain the smoothed mean of signal for
    % the importance model    
    par_KFS = IS_Model(par_NAIS, par_SV);
    [theta_smooth, ~, v, F] = KFS(y_star,par_KFS); % F is H


    %% LogLikelihood evaluation
    % loglikelihood of the approximation linear state space model via e.g. KF
    lng_y = -0.5*(n*log(2*pi) + sum(log(F)) + sum((v.^2)./F));
%   -0.5*(sum(log(F))) is gconst
%     lng_y = -0.5*(sum(log(F)) + sum((v.^2)./F)); <-- was
     
    %% Importance sampling:
    % simulation smoothing for linear state space models to sample S/2
    % trajectories for the signal e.g. via JSDK
   
    theta_sim = zeros(n,S);
    theta_sim(:,1:S/2) = IS_sim(S,n,theta_smooth,par_KFS, RND);
    % antithetic draws
    theta_sim(:,(S/2+1):S) = 2*repmat(theta_smooth,1,S/2) - theta_sim(:,1:S/2);
   
    % compute the logweights
    if (cont.err == 'n')
        lnP =  -0.5*(log(2*pi) + log(exp(theta_sim)) + repmat(y.^2,1,S)./exp(theta_sim));
    else % if (cont.err == 't')
        nu = par_SV(1,4);
        p_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2); %% pconst
        y2 = repmat(y.^2,1,S)./((nu-2).*exp(theta_sim));
        lnP = p_const - 0.5*(theta_sim + (nu+1)*log(1 + y2));
       % lnP = - 0.5*(theta_sim + (nu+1)*log(1 + y2));   %%   Y=logpdfkernel
    end
    lnG = repmat(a,1,S) + repmat(b,1,S).*theta_sim -0.5*repmat(C,1,S).*theta_sim.^2; %%  [X,gconst]=NAISfinalise

    logp = sum(lnP, 1);
    logg = sum(lnG, 1);
    lnw_s = logp - logg;
    max_lnw = max(lnw_s); % max log weight
    lnw_s = lnw_s - max_lnw; % robustification
    w_s = exp(lnw_s);
    lnw = mean(w_s);
    lnw = log(lnw);
    % compute the loglikelihood estimate
    lnL_hat =  lng_y + lnw + max_lnw;

%     lnL_hat = lnL_hat/n; % /n for stabilty?
%     fprintf('The estimated loglikelihood is: %6.4f.\n',lnL_hat)

    x = theta_sim(:,1);
    x = x';
    lnw_x = lnw_s(:,1) + max_lnw;

end