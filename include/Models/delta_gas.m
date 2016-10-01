function hessian = delta_gas(mu, hessian, T)
        % theta = [mu, omega, A, B]
        % constraints: 
        %               omega > 0; 
        %               0 < B < 1
        % [par_GAS_trans,~,~,~,~, hessian]= fminunc(posterior_gas_init);
 
%         if strcmp(cont.err,'n')
%             jaco_inv = diag([1, (1+exp(par_SV_trans(1,2)))*(1+exp(-par_SV_trans(1,2))), 1/(exp(par_SV_trans(1,3)))]);
%         else
%             jaco_inv = diag([1, (1+exp(par_SV_trans(1,2)))*(1+exp(-par_SV_trans(1,2))), 1/(exp(par_SV_trans(1,3))), 1/(exp(par_SV_trans(1,4)))]);
%         end
%         hessian_tr = jaco_inv*(n*hessian)*jaco_inv;
        omega = mu(1,2);
        B = mu(1,4);
        
        jaco_inv = diag([1, 1/(exp(omega)), 1, (1+exp(B))*(1+exp(-B))]);
        
        hessian = jaco_inv*(T*hessian)*jaco_inv;
end