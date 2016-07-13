function [loglik, f, param, score, scaledscore] = loglik_gas(theta, y, cont, GamMat)
%To enforce the volatility process to be positive, we model log volatility
% h_t = log(sigma2_t) rather than sigma2_t itself
    T = size(y,1);
    N = size(theta,1);
    
%     cont:
%     cont.err = 'n' / 't' // LATER ON: 'st' / 'mn' / 'mt' / 'ged' / 'egb2'
%     cont.link = 'id' / 'exp'
%     cont.scale = 'inv' / 'sqrtinv'

    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4); 
    if strcmp(cont.err,'t')
        nu = theta(:,5);
        rho = (nu-2)./nu;
    
        nu_con = (nu+1)./(nu-2);
        A = A.*((nu+3)./nu); 
        
        cond_ok = conditions_t_gas(omega, B, nu);
    elseif strcmp(cont.err,'n')
        cond_ok = conditions_t_gas(omega, B);
    end

    if strcmp(cont.link,'id') 
        fn_link = @(xx) xx;
    elseif strcmp(cont.link,'exp')
        fn_link = @(xx) exp(xx);        
    else
        Display('Incorrect link function.\n')
        Error('Incorrect link function.')
    end

    % initialize

    loglik = -Inf*ones(N,1);
    f = zeros(N,T);    
    param = zeros(N,T);
    
    for ii = 1:N
       
        pdf = zeros(T,1);
        
        if (cond_ok(ii,1)) % when all the parameter constraints are satisfied
            f(ii,1) = omega(ii,1)/(1-B(ii,1)); % unconditional variance to initialize h_1
            param(ii,1) = fn_link(f(ii,1));
            
            if  strcmp(cont.err,'t')
                pdf(1,1) = dmvt(y(1,1), mu(ii,1), rho(ii,1)*param(ii,1), nu(ii,1), GamMat);
                pdf(1,1) = log(pdf(1,1));
            elseif  strcmp(cont.err,'n')
                pdf(1,1) = -0.5*(log(2*pi) + log(param(ii,1)) + ((y(1,1)-mu(ii,1))^2)/param(ii,1));
            end
            
            
            
            for jj = 2:T
                C = 1 + ((y(jj-1,1)-mu(ii,1)).^2)/((nu(ii,1)-2)*f(jj-1,1));
                
                f(jj,1) = omega(ii,1) + A(ii,1)*(nu_con(ii,1)*((y(jj-1,1)-mu(ii,1)).^2)/C - f(jj-1,1)) ...
                            + B(ii,1)*f(jj-1,1);
                param = fn_link(f)        
                        
                if  strcmp(cont.err,'t')                  
                    pdf(jj,1) = dmvt(y(jj,1), mu(ii,1), rho(ii,1)*f(jj,1), nu(ii,1), GamMat);
                    pdf(jj,1) = log(pdf(jj,1));
                elseif  strcmp(cont.err,'n')
                    pdf(jj,1) = -0.5*(log(2*pi) + log(f(jj,1)) + ((y(jj,1)-mu(ii,1))^2)/f(jj,1));
                end
            end
            loglik(ii,1) = sum(pdf); 
        end
    end
       
    
    loglik = loglik/N;
end





function R = conditions_t_gas(omega, B, nu)
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c2 = ((B >= 0) & (B < 1));

    if (nargin == 2)
        R = (c1 & c2); 
    else
        c3 = (nu > 2);
        R = (c1 & c2 & c3); 
    end
end
