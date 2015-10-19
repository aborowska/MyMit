function d = posterior_Tgarch(theta, data, h1, L, hyper)
    % theta is Nxd, d=5, matrix of draws
    mu = theta(:,1);
    alpha0 = theta(:,2);
    alpha1 = theta(:,3);
    beta = theta(:,4);
    nu = theta(:,5);
    
    [N,~] = size(theta);
    T = length(data);
    %hyper = 0.01;
    prior = prior_Tgarch(mu, alpha0, alpha1, beta, nu, L, hyper);
    d = -Inf*ones(N,1); % density
    h = zeros(T,1); h(1,1) = h1;
    
    for ii = 1:N % row-wise through draws
        pdf = zeros(T,1);
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            rho = (nu(ii,1)-2)/nu(ii,1);
            for jj = 2:T
                h(jj,1) = alpha0(ii,1) + alpha1(ii,1)*(data(jj-1,1)-mu(ii,1))^2 + beta(ii,1)*h(jj-1,1);
                tmp = (data(jj,1)-mu(ii))/sqrt(rho*h(jj,1)); % standardized variable
                pdf(jj,1) = log(tpdf(tmp,nu(ii))); % standardized student-t
            end
            d(ii,1) = sum(pdf) + prior(ii,2); 
        end
    end
    if (~L)
        d = exp(d);
    end
end


function R = prior_Tgarch(mu, alpha0, alpha1, beta, nu, L, hyper)
    % non-informative (flat) priors for all paramters
    % a proper non-informative prior is used for nu-2 to avoid improper
    % posterior --> exponential prior with the restriction nu>2 to ensure that
    % the conditional variance is finite
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
%     c1 = ((alpha0 > 0) & (alpha0 < 1) & (alpha1 >= 0) & (alpha1 <= 1) &...
%         (beta >= 0) & (beta <= 1)); % to guarantee that the conditional variance will be positivr   
	c1 = ((alpha0 > 0) & (alpha0 < 1) & (alpha1 > 0) & (alpha1 < 1) &...
        (beta > 0) & (beta < 1)); % to guarantee that the conditional variance will be positivr    
    c2 = (nu > 2);
%     c3 = ((mu > -1) & (mu < 1));  %????
    
%     r1 = (c1 & c2 & c3);
    r1 = (c1 & c2);

    r2 = -Inf*ones(length(mu),1);
    % ????
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true) - 2); % exponential prior: p(nu)~exp(1) --> p(nu)=exp(-nu) from 2 to inf  
    if (~L)
        r2 = exp(r2);
    end
    R = [r1, r2];
end
