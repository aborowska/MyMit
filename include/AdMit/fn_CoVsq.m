function [CoVsq, d_CoVsq] =  fn_CoVsq(lnp, lnk, lnd)

     [Np, H] = size(lnk);
     p = exp(lnp);
     k = exp(lnk); 
     d = exp(lnd);  %     d = exp(lnd_up);
                % on pos = 1:H - evalutations of the draw from the first
                % compontent on the first, the second, ..., the new
                % component
                % ...
                % on pos = 1+(H-1)H : H^2 - evaluations of the draw from
                % the new component on the first, the second, ..., the new
                % component     
     tmp_eta = repmat(p, Np, H);
     tmp_k =  repmat(k, 1, H);
%      tmp_norm = sum(tmp_eta.*d,2);
     tmp_norm = tmp_eta.*d;
%      tmp_norm = repmat(tmp_norm, Np, 1);
     tmp_norm = tmp_norm * kron(eye(H),ones(H,1)); % the sum in the denominator of w(theta)
     
     sum_h = repmat(tmp_norm, 1, H);
     sum_d_nom = ((tmp_k.^2).*tmp_eta.*d./(sum_h.^3))* repmat(diag(ones(H,1)), H,1);
     sum_d_denom = (tmp_k.*tmp_eta.*d./(sum_h.^2))* repmat(diag(ones(H,1)), H,1);

     
     tmp_w = k./tmp_norm;
     tmp_eta = repmat(p, Np, 1);
     nom = sum(sum(tmp_eta.*tmp_w.^2))/Np;
     denom = sum(sum(tmp_eta.*tmp_w))/Np;
     
     d_nom = k.^2./(tmp_norm.^2) - 2.*sum_d_nom;
     d_nom = sum(d_nom,1)/Np;
     d_denom = k./tmp_norm - sum_d_denom;
     d_denom = 2*denom*sum(d_denom,1)/Np;
     
     CoVsq = -nom/(denom^2);
     d_CoVsq = -(d_nom.*denom - nom.*d_denom)./(denom^2);
     d_CoVsq = d_CoVsq';
end