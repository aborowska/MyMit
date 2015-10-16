function d = posterior_svt_whole(y, theta, x, prior_const, cont) 
% >>>>>>>>>>> ?NAIS?? 

    N = size(theta,1);
    
    nu = theta(:,4);
    eta = theta(:,5);
    eps = theta(:,6);   
   
    prior = prior_svt(theta(:,1:4),prior_const);

    par_SV = theta(:,1:4);

    if (N <= 2000)
        [~, lng_y, lnw_x] = NAIS_loglik(y, par_SV, par_NAIS, cont, x); 
    else 
        lng_y = zeros(N,1);
        lnw_x = zeros(N,1);
        for ii = 1:(N/1000)
            ind = (1:1000) + (ii-1)*1000; 
            par_NAIS_ind.b = par_NAIS.b(:,ind);
            par_NAIS_ind.C = par_NAIS.C(:,ind);
            [~, lng_y(ind,:), lnw_x(ind,:)] = NAIS_loglik(y, par_SV(ind,:), par_NAIS_ind, cont, x(ind,:)); 
        end
    end
   
    d = lng_y + lnw_x + prior;
    d = d + prior_const(1,1) - 0.5*eta.^2;
    d = d + prior_const(1,1) +  duvt(eps, nu, 1, true); 
end
