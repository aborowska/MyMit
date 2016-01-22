function theta_sim = IS_sim(y_star, v, F_inv, K, L, par_KFS)
% sample a signal trajectory via the simulation smoother
%     eps_sim = SimSmooth(v, F_inv, eps_smooth, K, L, par_KFS);
%     S = size(eps_sim,2);
%     theta_sim = kron(y_star,ones(1,S)) - eps_sim;
    [n, S] = size(y_star);
    RND = randn(n,S);
    theta_sim = SimSmooth(y_star, v, F_inv, K, L, par_KFS, RND);
%     S = size(eps_sim,2);
%     theta_sim = kron(y_star,ones(1,S)) - eps_sim;

end