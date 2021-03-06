function d = posterior_t_garch_hl(theta, data, S, VaR, L, hyper, GamMat)
    [N, hp] = size(theta);
    hp = hp - 4;
    alpha = theta(:,1);
    beta = theta(:,2);
    mu = theta(:,3);
    nu = theta(:,4);
    eps = theta(:,5:hp+4);
    theta = theta(:,1:4);
    
    T = size(data,1);
    y_T = data(T);
    ind = 2:T;

    % theta=[alpha, beta, mu, nu]
    h_T = volatility_t_garch(theta, data, S);
    y_hp = predict_t_garch(theta, y_T, S, h_T, hp, eps);
    PL_hp = fn_PL(y_hp);
    
    prior = prior_t_garch_hl(alpha, beta, mu, nu, PL_hp, VaR, L, hyper);
    d = -Inf*ones(N,1);
    h = zeros(T,1); 
    h(1,1) = S;
    omega = S*(1-alpha-beta); % variance targeting constraint
    rho = (nu-2)./nu;
    
    eps_pdf = duvt(eps, nu, hp, true); %log density
    eps_pdf = sum(eps_pdf, 2);
    
    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('posterior hl ii = %d\n',ii);
        end
        
        pdf = zeros(T,1);
        if (prior(ii,1)) % when all the parameter constraints are satisfied
%             h(1,1) = omega(ii,1)/(1-alpha(ii,1)-beta(ii,1)); % unconditional variance to initialize h_1
%             pdf(1,1) = dmvt(data(1,1), mu(ii,1), rho(ii,1)*h(1,1), nu(ii,1), GamMat);
%             pdf(1,1)= log(pdf(1,1));            
            
            h(ind) = omega(ii,1) + alpha(ii,1)*(data(ind-1,1)-mu(ii,1)).^2;
            for jj = ind
                h(jj,1) = h(jj,1) + beta(ii,1)*h(jj-1,1);
                pdf(jj,1) = dmvt(data(jj,1), mu(ii,1), rho(ii,1)*h(jj,1), nu(ii,1), GamMat);
                pdf(jj,1)= log(pdf(jj,1));
            end
            d(ii,1) = sum(pdf) +  prior(ii,2) + eps_pdf(ii,1); 
%             for hh = 1:hp
%                 d(ii,1) =  d(ii,1) + dmvt(eps(ii,hh), 0, 1, nu(ii,1));
%             end
        end
    end
    if (~L)
        d = exp(d);
    end
end


function R = prior_t_garch_hl(alpha, beta, mu, nu, PL_hp, VaR, L, hyper)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % uniform prior on mu on [-1,1]
    % exponential uninformative prior on nu with hyperparam hyper
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = ((mu > -1) & (mu < 1));
    c4 = (nu > 2);
    c5 = (PL_hp <= VaR);
    
    r1 = (c1 & c2 & c3 & c4 & c5);
    r2 = -Inf*ones(length(alpha),1);
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true) - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  
    if (~L)
        r2 = exp(r2);
    end
    R = [r1, r2];
end