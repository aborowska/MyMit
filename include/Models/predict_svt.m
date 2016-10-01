function [y_hp, eta_hp, eps_hp] = predict_svt(theta, x_T, H)
    [N, d] = size(theta);
    hp = (d - 4)/2;
   
    c = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    nu = theta(:,4);
    rho = (nu-2)./nu;
   
    if (hp == 0) % theta is only the model parameters
        eta_hp = randn(N,hp);               % prior: eta ~ N(0,1)
        eps_hp = trnd(repmat(nu,1,hp));     % prior: eps ~ t(0,1,nu)
    else % with given eta and eps      
        eta_hp = theta(:,5:4+hp);   
        eps_hp = theta(:,4+hp+1:d);   
        H = hp; % for robustness, adjust horizon length
    end    

    y_hp = zeros(N,H);    
     
    for jj = 1:H
        x_T = c + phi.*(x_T - c) + sqrt(sigma2).*eta_hp(:,jj);
        y_hp(:,jj) = sqrt(rho).*exp(0.5*x_T).*eps_hp(:,jj);
    end
end

