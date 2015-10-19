function R = prior_garch(theta, L)
    % uniform prior on 4-dimensional parameter vector theta=(mu,omega,alpah,beta)
    % [-1,1]X(0,1]X[0,1)X[0,1) with alpha+beta<1 for covariance stationarity
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    omega = theta(:,1);
    beta = theta(:,2);
    alpha = theta(:,3);
    mu = theta(:,4);
    c1 = ((omega > 0) & (omega < 1) & (beta >= 0) & (alpha >= 0));
    c2 = (beta + alpha < 1);
    c3 = ((mu > -1) & (mu < 1));
    r1 = (c1 & c2 & c3);
    r2 = -Inf*ones(length(omega),1);
    r2(r1==true) = 1;
    if (~L)
        r2 = exp(r2);
    end
    R = [r1, r2];
end
