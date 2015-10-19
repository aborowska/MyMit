function r2 = prior_sv_hl(data, theta, VaR, par_NAIS_init, prior_const, cont) 
    N = size(theta,1);
    T = size(data,1);
    S = cont.S;
    
    c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1);
    phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2);
    eps = theta(:,4);   % prior: eps ~ N(0,1)
    
    logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
    logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
    logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;

%     [~, ~, ~, b, C] = fn_nais_wgts(y, theta, cont, par_NAIS_init);
%     [lnL_hat, x, lnw_x, b, C] = fn_nais_wgts(y, theta, cont, par_NAIS_init);
    % sv one-step-ahead forecasts
    
    sigma2 = zeros(N,1);
    for ii = 1:N
        if (mod(ii,10) == 0)
            fprintf('\n',ii); 
            fprintf('nais_predict ii = %i\n',ii); 
            fprintf('\n',ii); 
        end
        RND = randn(T+1,S/2); 

        par_SV = theta(ii,1:3); 
        par_NAIS = NAIS_param(par_NAIS_init, data, par_SV, cont); % Algorithm 2: Efficient importance parameters via NAIS

%         par_NAIS.b = b(ii,:)';
%         par_NAIS.C = C(ii,:)';

        % sigma2=exp(x)=exp(c+alpha)
        [sigma2(ii,:), ~, ~] = NAIS_predict(par_SV, par_NAIS, data, cont, RND);
    end
    
    y = sqrt(sigma2).*eps;
    PL = fn_PL(y);
    c1 = (PL <= VaR);
 
    c2 = ((phi > 0) & (phi < 1));
    c3 = (s2 > 0);
       
    r1 = (c1 & c2 & c3);
    r2 = -Inf*ones(length(c),1);
    
    r2(r1==true) = logpdf_norm(c(r1==true));
    r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
    r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
    r2(r1==true) = r2(r1==true) + prior_const(1,1) - 0.5*(eps(r1==true)).^2;
    
%     if (~L)
%         r2 = exp(r2);
%     end
    
%     R = [r1, r2];
end