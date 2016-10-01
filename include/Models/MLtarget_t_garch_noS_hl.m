function [d, PL_hp] = MLtarget_t_garch_noS_hl(eps, theta_mle, h_mle, y_T, VaR)
    % theta is Nx5, matrix of draws
    [N, H] = size(eps);
  
    theta_mat = repmat(theta_mle,N,1);
    h_mat = repmat(h_mle,N,1);
    nu_mat = theta_mat(:,5);

    y_hp = predict_t_garch_noS(theta_mat, y_T, h_mat, H, eps);
    PL_hp = fn_PL(y_hp);

    prior = prior_t_garch_hl(N, PL_hp, VaR);     
    
    eps_pdf = duvt(eps, nu_mat, H, true); %log density
    eps_pdf = sum(eps_pdf, 2);
    
    d = prior(:,2) + eps_pdf; 

end


function R = prior_t_garch_hl(N, PL_hp, VaR)
    r1 = (PL_hp <= VaR);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = 0;
    
    R = [r1, r2];
end
