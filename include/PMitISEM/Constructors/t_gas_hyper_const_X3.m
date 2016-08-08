function [X, input_X] = t_gas_hyper_const_X3(input_X)
% theta: 1st col = alpha, 2+ cols = eps 
    theta = input_X.theta;    
    [N, h] = size(theta);
    H = h - 5;
    if (H > 0) % there is eps in the conditiponing set
        [y_H, ~, f_H] = predict_t_gas(theta(:,1:5), input_X.y_last, input_X.f_last, 1, theta(:,h));
        input_X.f_last = f_H;
        input_X.y_cum = input_X.y_cum + y_H;
        input_X.y_last = y_H;
        X = [ones(N,1), fn_PL(input_X.y_cum)];   
    else
        X = ones(N,1);
    end
end

