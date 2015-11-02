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