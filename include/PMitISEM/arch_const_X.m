function X = arch_const_X(theta, y_T, S)
% theta: 1st col = alpha, 2+ cols = eps 
    [N, h] = size(theta);
    if (h > 1) % there is eps
        
        hp = h-1;
        
%         y_hp = zeros(N,hp+1);    
%         y_hp(:,1) = y_T*ones(N,1); % the last observation in the first column
% 
%         f_ht = @(xx) sqrt(S + theta(:,1).*(xx.^2 - S)); % observation equation
% 
%         for ii = 2:(hp+1)
%             y_hp(:,ii) = eps(:,ii-1).*f_ht(y_hp(:,ii-1));
%         end        
%         X = [ones(N,1), fn_PL(y_hp(:,2:end))];

        y_hp = predict_arch(theta(:,1), y_T, S, hp, theta(:,2:h));
        X = [ones(N,1), fn_PL(y_hp)];       
    else
        X = ones(N,1);
    end
end

