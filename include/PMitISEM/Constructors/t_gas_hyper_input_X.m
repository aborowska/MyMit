function input_X = t_gas_hyper_input_X(theta, y)
    input_X.theta = theta;
%     input_X.h_T = volatility_t_garch_noS_mex(theta(:,1:5), data, S);
    input_X.f_T = volatility_t_gas_mex(theta(:,1:5), y);

end
% fn_input_X = @(xx) t_garch_input_X(xx, data, S);
% fn_input_X = @(xx) xx; % <-- WN and arch