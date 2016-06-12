function y_hp = predict_t_garch_new(theta, data, S, hp, eps)
    [N ,~] = size(theta);
    alpha = theta(:,1);
    beta = theta(:,2);
    mu = theta(:,3);
    nu = theta(:,4);
    
    y_T = data(end);
    omega = S*(1-alpha-beta); % variance targeting constraint
    rho = (nu-2)./nu;
    
    if (nargin == 4)
%         fprintf('hp = %i \n',hp);
        eps_hp = trnd(repmat(nu,1,hp));
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T*ones(N,1);

    h = zeros(N,hp+1); 
    h_T = volatility_t_garch_mex(theta, data, S);
    h(:,1) = h_T;
    
    for jj = 2:(hp+1)
         h(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);  
         y_hp(:,jj) = mu(:,1) + sqrt(rho(:,1).*h(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
end