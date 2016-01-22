function d = posterior_arch_hl(theta, data, S, VaR, L)
    [N,H] = size(theta);
    alpha = theta(:,1);
    eps = theta(:,2:H);
    H = H-1;
    T = length(data);

    y = data(T);
    ind = 2:T;

%     prior = prior_arch_hl(alpha, eps, y, S, VaR, L);
    prior = prior_arch_hl(alpha, eps, y, S, VaR); % in logs

    d = -Inf*ones(N,1);
    h = zeros(T,1); h(1,1) = S;
    omega = S*(1-alpha); % variance targeting constraint

    
    for ii = 1:N
        pdf = zeros(T,1);
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            h(ind) =  omega(ii,1) + alpha(ii,1)*(data(ind-1,1)).^2;
            pdf(ind) = -0.5*(log(2*pi) + log(h(ind)) + (data(ind).^2)./h(ind));
%             for jj = 2:T
%                 h(jj,1) = omega(ii,1) + alpha(ii,1)*(data(jj-1,1))^2;
%                 pdf(jj,1) = log(normpdf(data(jj,1),0,sqrt(h(jj,1))));
%             end
%             d(ii,1) = sum(pdf) + log(normpdf(eps(ii,1))) + prior(ii,2); 
            d(ii,1) = sum(pdf) + prior(ii,2)  - 0.5*(H*log(2*pi) + sum(eps(ii,1).^2,2)); 
        end
    end
    if (~L)
        d = exp(d - max(d));
    end
end

% function R = prior_arch_hl(alpha, eps, y, S, VaR, L)
function R = prior_arch_hl(alpha, eps, y, S, VaR)
    % uniform prior on a parameter alpha on [0,1)
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
%     c = 100*log(1+ VaR/100);
%     b = c*ones(length(alpha),1)./sqrt(S+(y^2-S).*alpha);
    H = size(eps,2) ;
    c1 = ((alpha >= 0) & (alpha < 1));
%     c2 = (eps <= b); 
    
    if (H == 1)
        h = S*(1-alpha) + (y^2)*alpha;
        c2 = (fn_PL(sqrt(h).*eps) <= VaR);
    else
        y_T1 = predict_arch(alpha, y, S, H, eps);
        c2 = (fn_PL(y_T1) <= VaR);
    end
    r1 = (c1 & c2);
    r2 = -Inf*ones(length(alpha),1);
    r2(r1==true) = 1;
%     if (~L)
%         r2 = exp(r2);
%     end
    R = [r1, r2];
end
