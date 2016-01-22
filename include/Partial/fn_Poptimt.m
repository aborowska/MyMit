function [mit_new, summary] = fn_Poptimt(theta, mit, w, cont, GamMat, X)
% optimizes the parameters of a H component mixture of k-variate student t
% densities using IS weighted EM algorithm, given N draws from the initial
% H component mixture
     
    % ISEM until convergence
    conv = 0; % convergence indicator
    iter = 0; % number of current iteratiosn
    loglik_old = - Inf; % previous loglik from ISEM
    mit_old = mit;
    H = length(mit_old.p);
        
    while (conv == 0) && (iter < cont.EM.maxit)
         shrinkmit = 1; % indicator for shringing mit density
        % apply/repeat optimization if mixture component is shrank
          while (shrinkmit)
            iter = iter + 1;
            fprintf('Iter in fn_Poptimt: %i\n',iter)
            % ISEM iteration
            mit_new = fn_PISEM(theta, mit_old, w, cont, GamMat, X);

            % shrink mit density only if H>1
             if (H == 1)
                 shrinkmit = 0;
    %            mit_old = mit_new; % ?????
             else
                 [shrinkmit, mit_old] = fn_shrink_mit(mit_new, cont.mit.tol_pr);
                 if (shrinkmit > 0) 
                    fprintf('Crash in fn.shrink.mit: %i \n',shrinkmit);
                 end
             end

             if (shrinkmit)
                 iter = 0;
             end
         end
         H = size(mit_new.p,2); % numebr of mixture components
         
         % new IS-weighted loglikelihood of the mixture
         tmp = dmvgt(theta, mit_new, true, GamMat);
         loglik_new = sum(w.*tmp);
         
         % convergence criteria
         if (iter > 1) && (abs((loglik_new - loglik_old)/loglik_old) <= cont.EM.tol)
             conv = 1;
         end
         loglik_old = loglik_new;
%          mit_old = mit_new; % without shrinking 

         if iter == cont.EM.maxit
            conv = 2; % maximum number of iterations reached
         end        
         summary = struct('EM_iter',iter,'EM_conv',conv);
    end
end 