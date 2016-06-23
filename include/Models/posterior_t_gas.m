function d = posterior_t_gas(theta, y, hyper, L, GamMat)
    % theta is Nx5, matrix of draws
    [N ,~] = size(theta);
    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4); 
    nu = theta(:,5);
    
    rho = (nu-2)./nu;
    
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);

    prior = prior_t_gas(N, omega, B, nu, hyper);
    
    T = size(y,1);
     
    d = -Inf*ones(N,1);
       
    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('posterior ii = %d\n',ii);
        end
        
        f = zeros(T,1);
        pdf = zeros(T,1);
        
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            f(1,1) = omega(ii,1)/(1-B(ii,1)); % unconditional variance to initialize h_1
            
            pdf(1,1) = dmvt(y(1,1), mu(ii,1), rho(ii,1)*f(1,1), nu(ii,1), GamMat);
            pdf(1,1) = log(pdf(1,1));
            
            for jj = 2:T
                C = 1 + ((y(jj-1,1)-mu(ii,1)).^2)/((nu(ii,1)-2)*f(jj-1,1));
                
                f(jj,1) = omega(ii,1) + A(ii,1)*(nu_con(ii,1)*((y(jj-1,1)-mu(ii,1)).^2)/C - f(jj-1,1)) ...
                            + B(ii,1)*f(jj-1,1);
                        
                pdf(jj,1) = dmvt(y(jj,1), mu(ii,1), rho(ii,1)*f(jj,1), nu(ii,1), GamMat);
                pdf(jj,1) = log(pdf(jj,1));
            end
            d(ii,1) = sum(pdf) + prior(ii,2); 
        end
    end
    
    if (~L)
        d = exp(d);
    end
end


function R = prior_t_gas(N, omega, B, nu, hyper)
    % uniform priors 
    
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    c4 = (nu > 2);
    
    r1 = (c1 & c3 & c4);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  

    R = [r1, r2];
end
