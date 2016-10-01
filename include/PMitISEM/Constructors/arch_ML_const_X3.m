function [X, input_X] = arch_ML_const_X3(input_X, S)
    theta = input_X.theta;
    [N, h] = size(theta);
    y_H = predict_arch(input_X.theta_mat, input_X.y_last, S, 1, theta(:,h));

    input_X.y_cum = input_X.y_cum + y_H;
    input_X.y_last = y_H;
    X = [ones(N,1), fn_PL(input_X.y_cum)];
end