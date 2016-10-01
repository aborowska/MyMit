function [X, input_X] = arch_const_X3(input_X, S)
% theta: 1st col = alpha, 2+ cols = eps 
    theta = input_X.theta;
    [N, h] = size(theta);
    if (h > 1) % there is eps
        hp = h-1;      
        y_H = predict_arch(theta(:,1), input_X.y_last, S, 1, theta(:,h));
        input_X.y_cum = input_X.y_cum + y_H;
        input_X.y_last = y_H;
        X = [ones(N,1), fn_PL(input_X.y_cum)];      
    else
        X = ones(N,1);
    end
end

