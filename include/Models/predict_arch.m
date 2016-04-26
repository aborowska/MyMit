function [y_hp, eps_hp] = predict_arch(alpha, y_T, S, hp, eps)
    % alpha - parameter vector
    % y_T - the last observation
    % S - sample variance
    % hp - horizion period
    % eps - (optional) error terms for prediction
    [N ,~] = size(alpha);
     
%     omega = S*(1-alpha); % variance targeting constraint
   
    if (nargin == 4)
        eps_hp = randn(N,hp);
    else % with given eps
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T*ones(N,1); % the last observation in the first column
    
    f_ht = @(xx) sqrt(S + alpha.*(xx.^2 - S)); % observation equation
    
    for ii = 2:(hp+1)
        y_hp(:,ii) = eps_hp(:,ii-1).*f_ht(y_hp(:,ii-1));
    end
%     h = zeros(N,hp+1); 
%     h(:,1) = omega;
    
%     y_hp_g = zeros(N,hp+1);    
%     y_hp_g(:,1) = y_T*ones(N,1);
    
%     h_g = zeros(N,hp+1); 
%     h_g(:,1) = omega;
    
%     ind = 2:hp+1;
    
%     h(:,ind) = repmat(omega(:,1),1,hp) + repmat(alpha(:,1),1,hp).*(y_hp(:,ind-1)).^2 ;  
%     y_hp(:,ind) = sqrt(h(:,ind)).*eps_hp;
%     
%     for jj = 2:(hp+1)
%         h(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)).^2;  
%         y_hp(:,jj) = sqrt(h(:,jj)).*eps_hp(:,jj-1);
%     end
 
    y_hp = y_hp(:,2:end);


end