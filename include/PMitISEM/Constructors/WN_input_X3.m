function input_X = WN_input_X3(theta)
    input_X.theta = theta;
    input_X.sigma = sqrt(theta(:,1));
    input_X.y_cum = zeros(size(theta,1),1);
end