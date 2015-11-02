function [y_hp, eps_hp, y_hp_g] = predict_arch(alpha, y_T, S, hp, eps)
    [N ,~] = size(alpha);
     
    omega = S*(1-alpha); % variance targeting constraint
   
    if (nargin == 4)
        eps_hp = randn(N,hp);
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T*ones(N,1);
    
    h = zeros(N,hp+1); 
    h(:,1) = omega;
    
    y_hp_g = zeros(N,hp+1);    
    y_hp_g(:,1) = y_T*ones(N,1);
    
    h_g = zeros(N,hp+1); 
    h_g(:,1) = omega;
    
    ind = 2:hp+1;
    
    h(:,ind) = repmat(omega(:,1),1,hp) + repmat(alpha(:,1),1,hp).*(y_hp(:,ind-1)).^2 ;  
    y_hp(:,ind) = sqrt(h(:,ind)).*eps_hp;
    
    for jj = 2:(hp+1)
        h_g(:,jj) = omega(:,1) + alpha(:,1).*(y_hp_g(:,jj-1)).^2;  
        y_hp_g(:,jj) = sqrt(h(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,ind);
    y_hp_g = y_hp_g(:,ind);

end