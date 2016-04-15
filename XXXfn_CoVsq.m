function CoVsq =  fn_CoVsq(lnp, lnkR, lndR)
%      [Np, H] = size(lnkR);
%      
%      p = exp(lnp);
%      k = exp(lnkR);
%      d = exp(lndR);
%                 % on pos = 1:H - evalutations of the draw from the first
%                 % compontent on the first, the second, ..., the new
%                 % component
%                 % ...
%                 % on pos = 1+(H-1)H : H^2 - evaluations of the draw from
%                 % the new component on the first, the second, ..., the new
%                 % component     
%      tmp_eta = repmat(p, Np, H);
% %      tmp_norm = sum(tmp_eta.*d,2);
%      tmp_norm = tmp_eta.*d;
% %      tmp_norm = repmat(tmp_norm, Np, 1);
%      tmp_norm = tmp_norm * kron(eye(H),ones(H,1)); % the sum in the denominator of w(theta)
% 
%      tmp_w = k./tmp_norm;
%      tmp_eta = repmat(p, Np, 1);
%      nom = sum(sum(tmp_eta.*tmp_w.^2))/Np;
%      denom = sum(sum(tmp_eta.*tmp_w))/Np;
% 
%      CoVsq = nom/(denom^2);
    [Np, H] = size(lnkR);
    
    k = lnkR';
    k = reshape(k, Np*H, 1);
    d = lndR';
    d = reshape(d, Np*H*H, 1);
    
    lnw = zeros(Np*H);
    nom = 0;
    denom = 0;
    m = 1; 
    n = 1;
    
    for ii = 1:Np
        for jj = 1:H
            tmp = 0;
            for kk = 1:H
                tmp = tmp + exp(lnp(kk) + d(m));                
                m = m + 1; % m - index of the log densitites
            end
            lnw(n) = k(n) - log(tmp); % log weight = log kernel/sum log densities
            nom = nom + exp(lnp(jj) + 2*lnw(n)); % probabilities jj range from 1 to H
            denom = denom + exp(lnp(jj) + lnw(n));
            n = n + 1; % n - index of the log kernels
        end
    end
    nom = nom/Np;
    denom = denom/Np;
    CoVsq = log(nom) - 2*log(denom);
end



% To minimize this function with the gradient provided, modify
% the function myfun so the gradient is the second output argument:
% function [f,g]= myfun(x)
% f = sin(x) + 3;
% g = cos(x);
% and indicate the gradient value is available by creating an options
% structure with OPTIONS.GradObj set to 'on' (using OPTIMSET):
% options = optimset('GradObj','on');
% x = fminunc(@myfun,4,options);