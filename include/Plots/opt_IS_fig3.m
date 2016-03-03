% optimal IS candidate standard normal profit/loss
figure(7)
xx = -5:0.01:5;
pl = normpdf(xx);
subplot(2,1,1)
plot(xx,pl)

p_bar = 0.01;  % 99% VaR
qq = norminv(p_bar); % hl subspace: xx<=qq
xx_inS = xx(xx<=qq);
xx_notS = xx(xx>qq);
pl_opt = [(1-p_bar)*normpdf(xx_inS), p_bar*normpdf(xx_notS)];
subplot(2,1,2)
plot(xx,pl_opt)

% xx_inS = linspace(-5,qq,500);
% xx_notS = linspace(qq,5,500);
% xx = [xx_inS, xx_notS];
% pl_opt = [(1-p_bar)*normpdf(xx_inS), p_bar*normpdf(xx_notS)];
% subplot(2,1,2)
% plot(xx,pl_opt)


figure(8)
xx = -5:0.01:5;
nu = 5;
pl = tpdf(xx,nu);
subplot(2,1,1)
plot(xx,pl)

p_bar = 0.01;  % 99% VaR
qq = tinv(p_bar,nu); % hl subspace: xx<=qq
xx_inS = xx(xx<=qq);
xx_notS = xx(xx>qq);
pl_opt = [(1-p_bar)*tpdf(xx_inS,nu), p_bar*tpdf(xx_notS,nu)];
subplot(2,1,2)
plot(xx,pl_opt)


% 
% yy = normpdf(xx)
% plot(sort(yy))
