function [theta, lnk, ind_red, x, lng_y, lnw_x, x_smooth] = fn_rmvgt_robust(N, mit, kernel, DUPA)
% robust sampling from mixture of multivariate t densities
% sampler are redrawn from mit if they correspond to a bad region with zero
% kernel density  (i.e. these with -Inf weights)        
        mu = mit.mu;
        Sigma = mit.Sigma;
        df = mit.df;
        p = mit.p;
        display('rmvgt2')

        x_smooth = [];
        
        if ((nargin == 4) && isa(DUPA,'double'))
            theta = DUPA; % theta=theta_init; cont=cont.nais; par_NAIS=par_NAIS_init;
        else
            theta = rmvgt2(N, mu, Sigma, df, p); 
        end
        
        % lnk - N vector of log-kernel evaluations at draws
        fprintf('\n'); 
        fprintf('kernel computation')
        fprintf('\n'); 
        
        if (nargin == 3)
            lnk = kernel(theta);
        else
%             [lnk, x, lng_y, lnw_x, x_smooth] =  kernel(theta);           
            [lnk, x, lng_y, lnw_x] =  kernel(theta);
        end
        
        ind_red = 0;
        while any(lnk == -Inf)
            ind_red = ind_red + 1;
            ind = find(lnk == -Inf);
            n_resamp = length(ind);
            fprintf('resampling %d draws.\n', n_resamp)
            draw_new = rmvgt2(n_resamp, mu, Sigma, df, p);
            theta(ind,:) = draw_new;

            if (nargin == 3)
                lnk(ind) = kernel(draw_new);
            else
%                 [lnk(ind), x(ind,:), lng_y(ind), lnw_x(ind), x_smooth(:,ind)] = kernel(draw_new);         
                [lnk(ind), x(ind,:), lng_y(ind), lnw_x(ind)] = kernel(draw_new);         
            end
        end
end 