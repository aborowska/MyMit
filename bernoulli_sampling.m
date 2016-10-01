function [S, F, theta] = bernoulli_sampling(T, K, theta_true, alpha_hyper, beta_hyper)
    if (nargin == 3)
        alpha_hyper = 1;
        beta_hyper = 1;         
    end
    S = zeros(1,K);
    F = zeros(1,K);
    theta = zeros(T,K);
    for t = 1:T
        A = S + alpha_hyper;
        B = F + beta_hyper;
        theta(t,:) = betarnd(A, B, 1, K);
        [~, arm] = max(theta(t,:));
        r = (rand < theta_true(1,arm));
        S(1,arm) = S(1,arm) + r;
        F(1,arm) = F(1,arm) + 1 - r;
    end
end
