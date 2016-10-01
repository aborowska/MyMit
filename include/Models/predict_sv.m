function [y_hp, eta_hp, eps_hp] = predict_sv(theta, x_T, H)
    [N, d] = size(theta);
    hp = (d - 3)/2;
    c = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);   
        
    if (hp == 0) % theta is only the model parameters
        eta_hp = randn(N,H);            % prior: eta ~ N(0,1)
        eps_hp = randn(N,H);            % prior: eps ~ N(0,1)
    else % with given eta and eps
        eta_hp = theta(:,4:3+hp);  
        eps_hp = theta(:,3+hp+1:d); 
        H = hp; % for robustness, adjust horizon length
    end

    y_hp = zeros(N,H);    
     
    for jj = 1:H
        x_T = c + phi.*(x_T - c) + sqrt(sigma2).*eta_hp(:,jj);
        y_hp(:,jj) = exp(0.5*x_T).*eps_hp(:,jj);
    end
end