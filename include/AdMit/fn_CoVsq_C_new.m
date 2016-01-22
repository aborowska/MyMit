function [CoVsq, d_CoVsq] =  fn_CoVsq_C_new(lnp, lnk, lnd)
     [Np, H] = size(lnk);
     % input lnp is the log of probabilities - we have unconstrained
     % optimisation (also negative values are allowed and not summing up to 1)
     % 1. transform back to have positivity and summability
     lnp = exp(lnp)/sum(exp(lnp));
     % 2. take log to have log-probabilities    
     p = log(lnp);        % p = log([0.9, 0.1]);

     k = lnk-repmat(median(lnk),Np,1);  % k = lnk_up(1:Np,:)-repmat(median(lnk_up(1:Np,:)),Np,1);    
     %lnk_up(1:Np,:), lnd_up
     d = lnd;  % d = lnd_up;   
                % on pos = 1:H - evalutations of the draw from the first
                % compontent on the first, the second, ..., the new
                % component
                % ...
                % on pos = 1+h*H : (h+1)*H - evaluation of the draw from
                % the h-th component on the first, the second, ..., the new
                % component
                % ...
                % on pos = 1+(H-1)H : H^2 - evaluations of the draw from
                % the new component on the first, the second, ..., the new
                % component    
     % replicate to have the same dimension as lnd           
     tmp_eta = repmat(p, Np, H);
     tmp_k =  repmat(k, 1, H);   
     % denomiator of w(theta), normalising constant for each draw (each
     % draw evaluated at all componets)
     tmp_norm = exp(tmp_eta + d);
     % sum evaluations of a draw at all components - dimention [NxH]
     tmp_norm = tmp_norm * kron(eye(H),ones(H,1)); % the sum in the denominator of w(theta)
     
     % replicate to have the same dimension as lnd           
     sum_h = repmat(tmp_norm, 1, H);
     
     sum_d_nom = exp(2*tmp_k + tmp_eta + d - 3*log(sum_h));
     % sum evaluations of all draws at a given component
     sum_d_nom  = sum_d_nom * repmat(diag(ones(H,1)), H,1);
     
     sum_d_denom = exp(tmp_k + tmp_eta + d - 2*log(sum_h));
     % sum evaluations of all draws at a given component
     sum_d_denom = sum_d_denom * repmat(diag(ones(H,1)), H,1);

     tmp_w = k - log(tmp_norm);
     tmp_eta = repmat(p, Np, 1);
     nom = exp(tmp_eta + 2*tmp_w);
     nom = sum(sum(nom))/Np;
     denom = exp(tmp_eta + tmp_w);
     denom = sum(sum(denom))/Np;
     
     d_nom = exp(2*k - 2*log(tmp_norm)) - 2.*sum_d_nom;
     d_nom = sum(d_nom,1)/Np;
     d_denom = exp(k - log(tmp_norm)) - sum_d_denom;
%      d_denom = 2*exp(log(denom) + log(sum(d_denom,1)))/Np;
     d_denom = 2*denom*sum(d_denom,1)/Np;
     
     
     CoVsq = -(log(nom) - 2*log(denom));
%      d_CoVsq = -(log(d_nom.*denom - nom.*d_denom) - 2*log(denom));
     d_CoVsq = -(d_nom.*denom - nom.*d_denom)./(denom^2);
     d_CoVsq = d_CoVsq';
end