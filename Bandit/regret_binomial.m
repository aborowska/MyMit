function regret = regret_binomial(theta)
    [G,K] =size(theta);
    theta_max = max(theta,[],2);
    ind = (theta == repmat(theta_max,1,K));
    w_a = sum(ind,1)/G;
    [~, opt_a] = max(w_a);
    theta_opt = theta(:,opt_a);
    regret = (theta_max - theta_opt)/theta_opt;
end