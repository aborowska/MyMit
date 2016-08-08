function [X, input_X] = WN_ML_const_X3(input_X)
    theta = input_X.theta;
    [N, h] = size(theta);

    y_H = input_X.sigma.*theta(:,h);
    input_X.y_cum = input_X.y_cum + y_H;
    X = [ones(N,1), fn_PL(input_X.y_cum)];    
end

