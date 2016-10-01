function input_X = t_gas_hyper_input_X3(theta, y)
    input_X.theta = theta;
%     input_X.h_T = volatility_t_garch_noS_mex(theta(:,1:5), data, S);
    input_X.f_last = volatility_t_gas_mex(theta(:,1:5), y);
    input_X.y_last = y(end)*ones(size(theta,1),1);
    input_X.y_cum = zeros(size(theta,1),1);    
end