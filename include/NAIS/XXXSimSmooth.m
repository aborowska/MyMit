function alpha_cond = SimSmooth(alpha_smooth_obs,par_KFS,RND)
% simulation smoothing for a linear Gaussian SSM
% i.e. conditional drawing given the observed sample y_obs
%     c = par_KFS.c; % mean adjustment
    H = par_KFS.H;
    Q = par_KFS.Q;
    T = par_KFS.T;
    R = par_KFS.R;
    Z = par_KFS.Z;
    
    n = length(alpha_smooth_obs);
    
    % Unconditional simulation
    eta_un = sqrt(Q).*RND(:,1);
    eps_un = sqrt(H).*RND(:,2);
    alpha_un = zeros(n+1,1);
    alpha_un(1,1) = alpha_smooth_obs(1,1);
    y_un = zeros(n,1);
    for ii = 1:n
        y_un(ii,1) = Z(ii,1)*alpha_un(ii,1) + eps_un(ii,1);
        alpha_un(ii+1,1) = T(ii,1)*alpha_un(ii,1) + R(ii,1)*eta_un(ii,1);
    end
    alpha_un = alpha_un(1:n,:);

    % Smoothed estimates
   alpha_smooth_un = KFS(y_un,par_KFS);

    % Conditional simulation
    alpha_cond = alpha_un + alpha_smooth_obs - alpha_smooth_un;
end