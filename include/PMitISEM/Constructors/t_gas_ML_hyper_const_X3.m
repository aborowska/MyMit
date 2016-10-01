function [X, input_X] = t_gas_ML_hyper_const_X3(input_X)
    theta = input_X.theta;
    [N, h] = size(theta);

    [y_H, ~, f_H] = predict_t_gas(input_X.theta_mat, input_X.y_last, input_X.f_last, 1, theta(:,h));
    input_X.f_last = f_H;
    input_X.y_cum = input_X.y_cum + y_H;
    input_X.y_last = y_H;
    X = [ones(N,1), fn_PL(input_X.y_cum)];   
end

