function [y_hp, eps_hp] = predict_t_garch(theta, y_T, S, h_T, hp, eps)
    [N ,~] = size(theta);
    alpha = theta(:,1);
    beta = theta(:,2);
    mu = theta(:,3);
    nu = theta(:,4);
      
    omega = S*(1-alpha-beta); % variance targeting constraint
    rho = (nu-2)./nu;
    
    if (nargin == 5)
         eps_hp = trnd(repmat(nu,1,hp));
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T*ones(N,1);

    h = zeros(N,hp+1); 
    h(:,1) = h_T;
    
    for jj = 2:(hp+1)
         h(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);  
         y_hp(:,jj) = mu(:,1) + sqrt(rho(:,1).*h(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
end