function [d, x, lng_y, lnw_x, x_smooth] = posterior_sv(y, theta, par_NAIS_init, prior_const, cont) 
    N = size(theta,1);
    T = size(y,1);
    x_smooth = [];
    prior = prior_sv(theta(:,1:3), prior_const);
  
    par_NAIS.b = -Inf*ones(T,N);
    par_NAIS.C = -Inf*ones(T,N);
%     x_smooth = zeros(T,N);
    for ii = 1:N
        if (mod(ii,100) == 0)
            fprintf('nais_param ii = %i\n',ii); 
        end
        par_SV = theta(ii,1:3); 
        % par_SV = [0.4982    0.9992    0.0002];
        % par_SV = [0.3630    0.9644    0.0562];

        % par_SV = mu_init;
        % par_SV = theta;
        % cont = cont.nais;
        % par_NAIS = par_NAIS_iter;
        % par_NAIS = par_NAIS_init;        
%         [par_NAIS_iter, x_smooth(:,ii)] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        [par_NAIS_iter] = NAIS_param(par_NAIS_init, y, par_SV, cont); % Efficient importance parameters via NAIS
        par_NAIS.b(:,ii) = par_NAIS_iter.b;
        par_NAIS.C(:,ii) = par_NAIS_iter.C;
    end
    par_SV = theta(:,1:3);
    fprintf('nais_loglik\n');
    [x, lng_y, lnw_x] = NAIS_loglik(y, par_SV, par_NAIS, cont); 

    d = lng_y + lnw_x + prior;
end

% function r2 = prior_sv(theta, prior_const)
%     c = theta(:,1);     % prior: c ~ normpdf(c, 0, 1);
%     phi = theta(:,2);   % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
%     s2 = theta(:,3);    % prior: 1/s2 ~ gampdf(1./s2, 5/2, 0.05/2);
% 
% logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
% logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
% logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) - 0.5*log(x) - 0.5*x;    
% 
%     c1 = ((phi > 0) & (phi < 1));
%     c2 = (s2 > 0);
%     r1 = (c1 & c2);
%     r2 = -Inf*ones(length(c),1);
%     
%     r2(r1==true) = logpdf_norm(c(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_beta((phi(r1==true)+1)/2);
% %     r2(r1==true) = r2(r1==true) + logpdf_gamma(1./s2(r1==true));
%     r2(r1==true) = r2(r1==true) + logpdf_invgamma(s2(r1==true));
% %     r2(r1==true) = r2(r1==true) + logpdf_chi2(s2(r1==true));
%     
% %     R = [r1, r2];
% end