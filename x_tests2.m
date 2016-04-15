theta2 = rmvgt2(1000, mit2.mu, mit2.Sigma, mit2.df, mit2.p);

% kernel = @(a) posterior_sv_x(y, a, par_NAIS_init, prior_const, cont.nais);
% [lnk, x, lng_y, lnw_x, eps_bar, eps_sim, C_sim, lnp_T] = kernel(theta2(ind-M/2,:));      


N = size(theta2,1);
T = size(y,1);
prior = prior_sv(theta2(:,1:3), prior_const);

par_NAIS.b = -Inf*ones(T,N);
par_NAIS.C = -Inf*ones(T,N);

theta_smooth = zeros(T,N);
for ii = 1:N
    if (mod(ii,100) == 0)
        fprintf('nais_param ii = %i\n',ii); 
    end
    par_SV = theta2(ii,1:3); 
    [par_NAIS_iter, theta_smooth(:,ii)] = NAIS_param(par_NAIS_init, y, par_SV, cont2.nais); % Efficient importance parameters via NAIS
    par_NAIS.b(:,ii) = par_NAIS_iter.b;
    par_NAIS.C(:,ii) = par_NAIS_iter.C;
end

par_SV = theta2;
theta_smooth = theta_smooth';

cont2.nais.HP = T-1-1;
fprintf('nais_loglik\n');
[x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] = NAIS_loglik_hl(y, par_SV, par_NAIS, cont2.nais); 

d = lng_y + lnw_x + prior ;

subplot(1,2,1)
hold all
for ii = 1:10
     plot(theta_smooth(ii,:))
end
hold off
title('theta smooth')
subplot(1,2,2)
hold all
for ii = 1:10
     plot(x(ii,:))
end
hold off
title('theta sim')

figure(21)
scatter(x(:,end-1), x(:,end))
x_end = x(:,end);
x_end = x_end(isfinite(x_end));
x_end1 = x(:,end-1);
x_end1 = x_end1(isfinite(x_end1));
cov(x_end,x_end1) % = [ 0.2596    0.2271;     0.2271    0.2406]
corr(x_end,x_end1) % = 0.9085


% CONSIDER HL SPACE
trans_prob = sum(lnp_T,2);
trans_prob_fin = trans_prob(isfinite(trans_prob));

trans_prob_T = lnp_T(:,end);
trans_prob_T_fin = trans_prob_T(isfinite(trans_prob_T));

lnw_x = lnw_x + sum(lnp_T,2);

c = theta2(:,1);
phi = theta2(:,2);
sigma2 = theta2(:,3);
eta = theta2(:,4);
eps = theta2(:,5);   
x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta;
y_h1 = exp(0.5*x_h1).*eps;    
PL = fn_PL(y_h1);
[PL_sort, ind_sort] = sort(PL);

    prior_hl = prior_sv_hl_in(theta, prior_const, PL, VaR);   
    ind = find(prior_hl ~= -Inf);
    d(ind) = lng_y(ind) + lnw_x(ind) + prior_hl(ind);

    
% CONSTRUCT ARTIFICIAL Xs    
mit_x.mu = 1.282975172394777;
mit_x.Sigma =  0.188166018837286;
mit_x.p = 1;
mit_x.df = 5;

x_T = rmvgt2(1000, mit_x.mu, mit_x.Sigma, mit_x.df, mit_x.p); 

x_const = c + phi.*(x_T - c) + sqrt(sigma2).*eta;
y_const= exp(0.5*x_const).*eps;    
PL_const = fn_PL(y_const);
[PL_const_sort, ind_const_sort] = sort(PL_const);

subplot(1,2,1)
plot(x_T) % sampled from mit_x
hold on
plot(x(:,end),'r') % from Sim Smoother
hold off
legend('from mit\_x','from SS')

subplot(1,2,2)
plot(PL_sort,'r')
hold on
plot(PL_const_sort,'b')
plot(VaR_prelim+0*PL_sort,'k')
hold off
legend('from SS','from mit\_x')

sum(PL_const_sort<VaR_prelim)
sum(PL_sort<VaR_prelim)










% What is the difference between the NAIS parameters for a shorter y??
par_NAIS_init_short.b = zeros(T-10,1);
par_NAIS_init_short.C = ones(T-10,1);

par_NAIS_short.b = -Inf*ones(T-10,N);
par_NAIS_short.C = -Inf*ones(T-10,N);
theta_smooth_short = zeros(T-10,N);

for ii = 1:N
    if (mod(ii,100) == 0)
        fprintf('nais_param ii = %i\n',ii); 
    end
    par_SV = theta2(ii,1:3); 
    [par_NAIS_iter, theta_smooth_short(:,ii)] = NAIS_param(par_NAIS_init_short, y(1:T-10), par_SV, cont2.nais); % Efficient importance parameters via NAIS
    par_NAIS_short.b(:,ii) = par_NAIS_iter.b;
    par_NAIS_short.C(:,ii) = par_NAIS_iter.C;
end

nais_diff_b = abs(par_NAIS_short.b - par_NAIS.b(1:T-10,:))';
nais_diff_C = abs(par_NAIS_short.C - par_NAIS.C(1:T-10,:))';

Mb = max(nais_diff_b);
MC = max(nais_diff_C);

Ab = mean(nais_diff_b);
AC = mean(nais_diff_C);

subplot(3,1,1)
plot(1:(T-10),y(1:T-10))
title('y ')

subplot(3,1,2)
hold on
plot(1:(T-10),Mb)
% title('max abs diff b')
plot(1:(T-10),Ab,'r')
hold off
title('diff b')
legend('max','avg')

subplot(3,1,3)
hold on
plot(1:(T-10),MC)
% title('max abs diff C')
plot(1:(T-10),AC,'r')
hold off
title('diff C')
legend('max','avg')

comment = {'Max and avg abs difference','between NAIS parameters','comuted for the whole y and','without the last 10 observation'};
% comment = {'Max abs difference','between NAIS parameters','comuted for the whole y and','without the last 10 observation'};
dim = [0.65  0.9  0.1 0.1];
% dim - four-element vector of the form [x y w h]. 
% x, y - the position; w, h - the size.
annotation('textbox',dim,'String',comment,'BackgroundColor','y');%,'FitBoxToText','on');





