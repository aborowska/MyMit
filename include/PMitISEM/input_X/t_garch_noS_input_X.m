function input_X = t_garch_noS_input_X(theta, data, S)
    input_X.theta = theta;
    input_X.h_T = volatility_t_garch_noS_mex(theta(:,1:5), data, S);
end
% fn_input_X = @(xx) t_garch_input_X(xx, data, S);
% fn_input_X = @(xx) xx; % <-- WN and arch