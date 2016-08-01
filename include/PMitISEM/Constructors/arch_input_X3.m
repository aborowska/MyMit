function input_X = arch_input_X3(theta, y)
    input_X.theta = theta;
    input_X.y_last = y(end)*ones(size(theta,1),1);
    input_X.y_cum = zeros(size(theta,1),1);    
end