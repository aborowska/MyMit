% SSM and VaR
% http://homepages.gold.ac.uk/nikolaev/varest.htm

% what the last 10 states look like in direct approach
% in particular: those leading to high losses
cont.nais.HP = 9;

theta1 = rmvgt2(M, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
eta_h1_1 = randn(M,1);
eps_h1_1 = randn(M,1);

theta1 = [theta1, eta_h1_1, eps_h1_1];

kernel = @(a) posterior_sv_x(y, a, par_NAIS_init, prior_const, cont.nais);

lnk = zeros(M,1);
x = zeros(M,cont.nais.HP+1);
lng_y = zeros(M,1);
lnw_x = zeros(M,1);
eps_bar = zeros(M,1);
eps_sim = zeros(M,1);
C_sim = zeros(M,1);
lnp_T = zeros(M,1);

for ii = 1:(M/1000)
    fprintf('ii = %i\n',ii)
    ind = (1:1000) + (ii-1)*1000; 
    [lnk(ind,:), x(ind,:), lng_y(ind,:), lnw_x(ind,:), eps_bar(ind,:), eps_sim(ind,:), C_sim(ind,:), lnp_T(ind,:)] = kernel(theta1(ind,:));      
end

theta_opt = [theta1, x(:,end)]; % extended draw - includes the last state

lnk = lnk + 2*prior_const(1,1) - 0.5*(theta_opt(:,4)).^2 - 0.5*(theta_opt(:,5)).^2;
exp_lnd = exp(2*prior_const(1,1) - 0.5*(theta_opt(:,4)).^2  - 0.5*(theta_opt(:,5)).^2 + dmvgt(theta_opt(:,1:3), mit1, true, GamMat));
lnd_opt = log(exp_lnd);

w_opt = fn_ISwgts(lnk, lnd_opt, false);

c_opt = theta_opt(:,1);
phi_opt = theta_opt(:,2);
sigma2_opt = theta_opt(:,3);

eta_opt = theta_opt(:,4);
eps_opt = theta_opt(:,5);  
x_opt_h1 = c_opt + phi_opt.*(x(:,end) - c_opt) + sqrt(sigma2_opt).*eta_opt;
y_opt_h1 = exp(0.5*x_opt_h1).*eps_opt;


dens = struct('y',y_opt_h1,'w',w_opt,'p_bar',p_bar);
IS_estim = fn_PL(dens, 1);
VaR_IS = IS_estim(1,1);
ES_IS = IS_estim(1,2);

save(['results/last10.mat'], 'mit1', 'theta_opt', 'x', 'lng_y','lnw_x',...
'eps_bar','eps_sim','C_sim','lnp_T', 'lnk', 'lnp_T', 'lnd_opt', 'w_opt', 'y_opt_h1', ...
 'cont',  'p_bar', 'M', 'VaR_prelim','VaR_IS', 'ES_IS');


ind_y = (imag(y_opt_h1)==0); 
y_opt_h1 = y_opt_h1(ind_y);
w_opt = w_opt(ind_y);
lnp_T = lnp_T(ind_y);
x = x(ind_y,:);

f_pl = @(aa) 100*(exp(aa/100) - 1); 
PL_opt_h1 = f_pl(sum(y_opt_h1,2));
[PL_opt_h1, ind] = sort(PL_opt_h1);
lnp_T = lnp_T(ind);
x = x(ind,:);

% take below var prelim or is
ind2 = (PL_opt_h1 <= VaR_IS); 
M2 = sum(ind2);
x2 = x(ind2,:);


ind1 = ((PL_opt_h1 > VaR_IS) & (PL_opt_h1<0)); % losses above var
M1 = sum(ind1);
ind11 = (ceil(M1/2) - floor(M2/2)):(ceil(M1/2) + floor(M2/2)); % 'moderate' losses
ind111 = zeros(length(ind2),1);
ind111(ind11,1) = 1;
ind111 = logical(ind111);
x1 = x(ind111,:);


sum_x1 = sum(x1,2);
sum_x2 = sum(x2,2);


set(gcf,'units','normalized','outerposition',[0 0 0.25 1]);   
set(gcf,'defaulttextinterpreter','latex');

subplot(3,1,1)
hist([sum_x1,sum_x2])
legend('moderate losses','high losses')
title('Sum of the last 10 x`es (log-volatilities), samples of 105 elements')    
plotTickLatex2D;

subplot(3,1,2)
hist([sum(x1(:,(end-2):end),2),sum(x2(:,(end-2):end),2)])
legend('moderate losses','high losses')
title('Sum of the last 3 x`es (log-volatilities), samples of 105 elements')    
plotTickLatex2D;

subplot(3,1,3)
hist([x1(:,end),x2(:,end)])
legend('moderate losses','high losses')
title('The last x (today`s log-volatility), samples of 105 elements')    
plotTickLatex2D;

figure(2)
hold on
plot(PL_opt_h1,'k')
plot(PL_opt_h1(ind2),'r','LineWidth',2)
plot(PL_opt_h1(ind111),'b','LineWidth',2)
hold off
