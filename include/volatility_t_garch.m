function h_T = volatility_t_garch(theta, data, S)
% function h_T = volatility_t_garch(theta, data, S)

    [N ,~] = size(theta);
    alpha = theta(:,1);
    beta = theta(:,2);
    mu = theta(:,3);
      
    omega = S*(1-alpha-beta); % variance targeting constraint

    T = size(data,1);
    ind = 2:T;
    data = data';
    
    h = zeros(N,T); 
    h(:,1) = S*ones(N,1);
    
    h(:,ind) = repmat(omega(:,1),1,T-1) + repmat(alpha(:,1),1,T-1).*(repmat(data(1,ind-1),N,1)-repmat(mu(:,1),1,T-1)).^2;
    for jj = ind
        h(:,jj) = h(:,jj)  + beta(:,1).*h(:,jj-1) ;
    end
    h_T = h(:,T);
end