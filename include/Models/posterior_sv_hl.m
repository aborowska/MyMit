function [lnk, x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] = posterior_sv_hl(y, theta, VaR, par_NAIS_init, prior_const, cont) 
    [N, d] = size(theta);
    hp = (d - 3)/2;
    T = size(y,1);
    
    b = -Inf*ones(T,N);
    C = -Inf*ones(T,N);
  
    lnk = -Inf*ones(N,1);

    for ii = 1:N      
        if (mod(ii,100) == 0)
            fprintf('Posterior sv_hl: nais_loglik ii = %i\n',ii); 
        end
        par_SV = theta(ii,1:3); 
        [par_NAIS_iter] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        b(:,ii) = par_NAIS_iter.b;
        C(:,ii) = par_NAIS_iter.C;
    end
    par_NAIS.b = b;
    par_NAIS.C = C;
    par_SV = theta(:,1:3);
    [x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] = NAIS_loglik(y, par_SV, par_NAIS, cont); 

    y_hp = predict_sv(theta, x(:,end), hp);
    PL = fn_PL(y_hp);

    prior_hl = prior_sv_hl(theta, prior_const, PL, VaR);   

    ind = find(prior_hl ~= -Inf);
    lnk(ind) = lng_y(ind) + lnw_x(ind) + prior_hl(ind);
end

function r2 = prior_sv_hl(theta, prior_const, PL, VaR)
    d = size(theta,2);
    
    c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1);
    phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2);

    if (d > 3)
        hp = (d - 3)/2;        
        eta = theta(:,4:3+hp);   % prior: eta ~ N(0,1)
        eps = theta(:,3+hp+1:d);   % prior: eps ~ N(0,1)
    end    
    
    logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
    logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
    % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
    logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
    % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) - 0.5*log(x) - 0.5*x;    

    c1 = (PL <= VaR);
    c2 = ((phi > 0) & (phi < 1));
    c3 = (s2 > 5e-3);   
  
    r1 = (c1 & c2 & c3);
    r2 = -Inf*ones(length(eta),1);
    
    r2(r1==true) = logpdf_norm(c(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
%     r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_invgamma(s2(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));

    if (d > 3)
        r2(r1==true) = r2(r1==true) + hp*prior_const(1,1) - 0.5*sum((eta(r1==true)).^2,2);
        r2(r1==true) = r2(r1==true) + hp*prior_const(1,1) - 0.5*sum((eps(r1==true)).^2,2);
    end
end