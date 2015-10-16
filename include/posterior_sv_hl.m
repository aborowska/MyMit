function [d, x, lng_y, lnw_x, x_smooth] = posterior_sv_hl(y, theta, VaR, par_NAIS_init, prior_const, cont) 
    N = size(theta,1);
    T = size(y,1);
    
    x_smooth = [];
    
    c = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
%     nu = theta(:,4);
    eta = theta(:,4);
    eps = theta(:,5);   
    
    b = -Inf*ones(T,N);
    C = -Inf*ones(T,N);
  
    d = -Inf*ones(N,1);
%     x_smooth = zeros(T,N);
    for ii = 1:N      
        if (mod(ii,100) == 0)
            fprintf('\n'); 
            fprintf('nais_loglik ii = %i\n',ii); 
            fprintf('\n'); 
        end
        par_SV = theta(ii,1:3); 
%         [par_NAIS_iter, x_smooth(:,ii)] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        [par_NAIS_iter] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        b(:,ii) = par_NAIS_iter.b;
        C(:,ii) = par_NAIS_iter.C;
    end
    par_NAIS.b = b;
    par_NAIS.C = C;
    par_SV = theta(:,1:3);
    [x, lng_y, lnw_x] = NAIS_loglik(y, par_SV, par_NAIS, cont); 

    x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta;
    y_h1 = exp(0.5*x_h1).*eps;    
    PL = fn_PL(y_h1);

    prior_hl = prior_sv_hl_in(theta, prior_const, PL, VaR);   

    ind = find(prior_hl ~= -Inf);
    d(ind) = lng_y(ind) + lnw_x(ind) + prior_hl(ind);

    par_NAIS.b = -Inf*ones(T,N);
    par_NAIS.C = -Inf*ones(T,N);
    par_NAIS.b(:,ind) = b(:,ind);
    par_NAIS.C(:,ind) = C(:,ind);

end

function r2 = prior_sv_hl_in(theta, prior_const, PL, VaR)
    c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1);
    phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2);
    eta = theta(:,4);   % prior: eta ~ N(0,1)
    eps = theta(:,5);   % prior: eps ~ N(0,1)
    
    logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
    logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
    % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
    logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
    % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) - 0.5*log(x) - 0.5*x;    

    c1 = (PL <= VaR);
    c2 = ((phi > 0) & (phi < 1));
%     c3 = (s2 > 0);
    c3 = (s2 > 5e-3);   
  
    r1 = (c1 & c2 & c3);
    r2 = -Inf*ones(length(eta),1);
    
    r2(r1==true) = logpdf_norm(c(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
%     r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_invgamma(s2(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));
    r2(r1==true) = r2(r1==true) + prior_const(1,1) - 0.5*(eta(r1==true)).^2;
    r2(r1==true) = r2(r1==true) + prior_const(1,1) - 0.5*(eps(r1==true)).^2;

%     R = [r1, r2];
end