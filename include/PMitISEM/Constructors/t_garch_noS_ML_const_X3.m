function [X, input_X] = t_garch_noS_ML_const_X3(input_X)
    theta = input_X.theta;
    [N, h] = size(theta);

    [y_H, ~, h_H] = predict_t_garch_noS(input_X.theta_mat, input_X.y_last, input_X.h_last, 1, theta(:,h));  
    input_X.h_last = h_H;
    input_X.y_cum = input_X.y_cum + y_H;
    input_X.y_last = y_H;
    X = [ones(N,1), fn_PL(input_X.y_cum)];   
end

