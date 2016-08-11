function input_X = t_garch_noS_ML_input_X3(theta, theta_mle, h_mle, y_T)
    N = size(theta,1);
    input_X.theta = theta;
    input_X.theta_mat = repmat(theta_mle,N,1);
    input_X.h_last = h_mle*ones(N,1);
    input_X.y_last = y_T*ones(N,1);
    input_X.y_cum = zeros(N,1);    
end