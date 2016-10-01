function [X, input_X] = WN_const_X3(input_X)
% theta: 1st col = alpha, 2+ cols = eps 
    theta = input_X.theta;
    [N, h] = size(theta);
    if (h > 1) % there is eps      
%         y_H = sqrt(theta(:,1)).*theta(:,h);
        y_H = input_X.sigma.*theta(:,h);
        input_X.y_cum = input_X.y_cum + y_H;
        X = [ones(N,1), fn_PL(input_X.y_cum)];      
    else
        X = ones(N,1);
    end
end

