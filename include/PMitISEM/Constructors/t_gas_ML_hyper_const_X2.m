function X = t_gas_ML_hyper_const_X2(eps, y_T, theta_mle, f_mle)
    [N, H] = size(eps);
    
%     if (H > 1) % there is eps in the conditiponing set
        theta_mat = repmat(theta_mle,N,1);
        f_mat = repmat(f_mle,N,1);
        y_H = predict_t_gas(theta_mat, y_T, f_mat, H, eps);
        X = [ones(N,1), fn_PL(y_H)];   
%     else
%         X = ones(N,1);
%     end
end

