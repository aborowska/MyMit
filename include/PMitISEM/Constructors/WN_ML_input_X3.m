function input_X = WN_ML_input_X3(theta, sigma2_used)
    input_X.theta = theta;
    N = size(theta,1);
    input_X.sigma = sqrt(sigma2_used)*ones(N,1);
    input_X.y_cum = zeros(N,1);
end