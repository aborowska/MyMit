function X = arch_const_X(theta, y_T, S)
% theta: 1st col = alpha, 2+ cols = eps 
    [N, h] = size(theta);
    if (h > 1) % there is eps
        hp = h-1;      
        y_hp = predict_arch(theta(:,1), y_T, S, hp, theta(:,2:h));
        X = [ones(N,1), fn_PL(y_hp)];       
    else
        X = ones(N,1);
    end
end

