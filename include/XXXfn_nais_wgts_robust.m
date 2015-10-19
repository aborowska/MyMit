function [theta, prior, lnd, lnL_hat, x, lnw_x, ind_red, n_resampl_total] = fn_nais_wgts_robust(y, theta, prior, lnd, lnL_hat, x, lnw_x, mit, kernel, cont, GamMat, par_NAIS_init)    
    ind_red = 0;
    n_resampl_total = 0;
    display('resampling')
    
    while any(~isfinite(lnL_hat)) %any(lnL_hat == -Inf)
        ind_red = ind_red + 1;
        ind = find(~isfinite(lnL_hat)) ; %(lnk == -Inf);
        n_resampl = length(ind);
        n_resampl_total = n_resampl_total + n_resampl;

        fprintf('resampling ind_red = %i, n_resmapl = %i\n',ind_red, n_resampl); 
        for ii = ind
            [theta_resampl, prior_resampl, ~] = fn_rmvgt_robust(n_resampl, mit, kernel);
            lnd_resampl = dmvgt(theta_resampl, mit, true, GamMat);
            [lnL_hat_resampl, x_resampl, lnw_x_resampl] = fn_nais_wgts(y, theta_resampl, cont, par_NAIS_init);
        end
        theta(ind,:) = theta_resampl;
        prior(ind,:) = prior_resampl;
        lnd(ind,:) = lnd_resampl;
        lnL_hat(ind,:) = lnL_hat_resampl;
        x(ind,:) = x_resampl;
        lnw_x(ind,:) = lnw_x_resampl;               
    end
end