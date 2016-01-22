function  [d, x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] = posterior_sv(y, theta, par_NAIS_init, prior_const, cont) 
    N = size(theta,1);
    T = size(y,1);
    prior = prior_sv(theta(:,1:3), prior_const);
  
    par_NAIS.b = -Inf*ones(T,N);
    par_NAIS.C = -Inf*ones(T,N);
    for ii = 1:N
        if (mod(ii,100) == 0)
            fprintf('nais_param ii = %i\n',ii); 
        end
        par_SV = theta(ii,1:3); 
        [par_NAIS_iter] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        par_NAIS.b(:,ii) = par_NAIS_iter.b;
        par_NAIS.C(:,ii) = par_NAIS_iter.C;

    end
    par_SV = theta(:,1:3);
    fprintf('nais_loglik\n');
    [x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] = NAIS_loglik(y, par_SV, par_NAIS, cont); 

    d = lng_y + lnw_x + prior;
end