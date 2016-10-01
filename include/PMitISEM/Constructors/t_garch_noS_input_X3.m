function input_X = t_garch_noS_input_X3(theta, data, S)
    input_X.theta = theta;
    input_X.h_last = volatility_t_garch_noS_mex(theta(:,1:5), data, S);
    input_X.y_last = data(end)*ones(size(theta,1),1);
    input_X.y_cum = zeros(size(theta,1),1);
end