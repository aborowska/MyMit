function X = t_gas_hyper_const_X2(input_X, y)
% theta: 1st col = alpha, 2+ cols = eps 
    if ~isstruct(input_X)
        theta = input_X;
    else
        theta = input_X.theta;
    end    
     
    [N, h] = size(theta);
    H = h - 5;
    if (H > 0) % there is eps in the conditiponing set
        y_T = y(end);
        if ~isstruct(input_X)
            f_T = volatility_t_gas_mex(theta(:,1:5), y);
%             h_T = volatility_t_garch_noS_mex(theta(:,1:5), data, S);
        else
            f_T = input_X.f_T;
        end 
%         y_H = predict_t_garch_noS(theta(:,1:5), y_T, S, h_T, H, theta(:,6:h));
        y_H = predict_t_gas(theta(:,1:5), y_T, f_T, H, theta(:,6:h));
        X = [ones(N,1), fn_PL(y_H)];   
    else
        X = ones(N,1);
    end
end

