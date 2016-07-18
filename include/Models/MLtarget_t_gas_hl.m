function [d, PL_hp] = MLtarget_t_gas_hl(eps, theta_mle, f_mle, y_T, VaR)
    % theta is Nx5, matrix of draws
    [N, H] = size(eps);
  
    theta_mat = repmat(theta_mle,N,1);
    f_mat = repmat(f_mle,N,1);
    nu_mat = theta_mat(:,5);
    
    y_hp = predict_t_gas(theta_mat, y_T, f_mat, H, eps);
    PL_hp = fn_PL(y_hp);

    prior = prior_t_gas_hl(N, PL_hp, VaR);     
    
    eps_pdf = duvt(eps, nu_mat, H, true); %log density
    eps_pdf = sum(eps_pdf, 2);
    
    d = prior(:,2) + eps_pdf; 

end


function R = prior_t_gas_hl(N, PL_hp, VaR)
    r1 = (PL_hp <= VaR);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = 0;
    
    R = [r1, r2];
end
