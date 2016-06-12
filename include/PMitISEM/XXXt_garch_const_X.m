function X = t_garch_const_X(theta, data, S)
% theta: 1st col = alpha, 2+ cols = eps 
    [N, h] = size(theta);
    H = h - 4;
    if (H > 0) % there is eps in the conditiponing set
        y_T = data(end);
%         filename = ['temp_h_T_',num2str(N),'.mat'];
%         if (H == 1) % the last volatility h_T to be computed 
            h_T = volatility_t_garch_mex(theta(:,1:4), data, S);
            % save h_T
%             save(filename,'-v7','h_T');
            y_H = predict_t_garch(theta(:,1:4), y_T, S, h_T, H, theta(:,5:h));
%         else % h_T has been computed
%             % load h_T
%             if (exist(filename, 'file') == 2)
%                 load(filename,'h_T');
%             else
%                  h_T = volatility_t_garch_mex(theta(:,1:4), data, S);
%             end
%             y_H = predict_t_garch(theta, y_T, S, h_T, H, theta(:,5:h));
%         end
        X = [ones(N,1), fn_PL(y_H)];   
    else
        X = ones(N,1);
    end
end

