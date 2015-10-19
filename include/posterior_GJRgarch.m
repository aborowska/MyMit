function d = posterior_GJRgarch(theta, data, S, L, hyper, GamMat)
    % theta is Nx4, matrix of draws
    [N ,~] = size(theta);
    alpha = theta(:,1);
    beta = theta(:,2);
    mu = theta(:,3);
    nu = theta(:,4);
    gamma = theta(:,5);

    prior = prior_GJRgarch(N, alpha, beta, mu, nu, gamma, L, hyper);    T = size(data,1);
    ind = 1:T-1;
    
    d = -Inf*ones(N,1);
    h = zeros(T,1); h(1,1) = S;
    omega = S*(1-alpha-beta); % variance targeting constraint
    rho = (nu-2)./nu;
       
    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('posterior ii = %d\n',ii);
        end
        
        pdf = zeros(T,1);
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            if (data(ind,1)>mu(ii,1))
                h(2:T) = omega(ii,1) + alpha(ii,1)*(data(ind,1)-mu(ii,1)).^2;
            else
                h(2:T) = omega(ii,1) + (alpha(ii,1)+gamma(ii,1))*(data(ind,1)-mu(ii,1)).^2;
            end
            
            for jj = 2:T
                h(jj,1) = h(jj,1) + beta(ii,1)*h(jj-1,1);
                pdf(jj,1) = dmvt(data(jj,1), mu(ii,1), rho(ii,1)*h(jj,1), nu(ii,1), GamMat);
                pdf(jj,1)= log(pdf(jj,1));
            end
            d(ii,1) = sum(pdf) + prior(ii,2); 
        end
    end
    if (~L)
        d = exp(d);
    end
end


function R = prior_GJRgarch(N, alpha, beta, mu, nu, gamma, L, hyper)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % uniform prior on mu on [-1,1]
    % exponential uninformative prior on nu with hyperparam hyper
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point

    c1 = ((alpha >= 0) & (beta >= 0) & (gamma >= 0));
    c2 = ((alpha + gamma)/2 + beta < 1);
    c3 = ((alpha + gamma)/2 > 0);
    c4 = (nu > 2);
    
    r1 = (c1 & c2 & c3 & c4);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  
    if (~L)
        r2 = exp(r2);
    end
    R = [r1, r2];
end
