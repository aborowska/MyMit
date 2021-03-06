function X = t_garch_const_X2(input_X, data, S)
% theta: 1st col = alpha, 2+ cols = eps 
    if ~isstruct(input_X)
        theta = input_X;
    else
        theta = input_X.theta;
    end    
     
    [N, h] = size(theta);
    H = h - 4;
    if (H > 0) % there is eps in the conditiponing set
        y_T = data(end);
        if ~isstruct(input_X)
            h_T = volatility_t_garch_mex(theta(:,1:4), data, S);
        else
            h_T = input_X.h_T;
        end 
        y_H = predict_t_garch(theta(:,1:4), y_T, S, h_T, H, theta(:,5:h));
        X = [ones(N,1), fn_PL(y_H)];   
    else
        X = ones(N,1);
    end
end

