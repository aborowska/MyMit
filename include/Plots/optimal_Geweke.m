x = -5:0.01:5;
y = tpdf(x,5);
VaR = tinv(0.01,5);


y_hl = y(x<=VaR);
y_hl = y_hl/sum(y_hl);
y_out = y(x>VaR);
y_out = y_out/sum(y_out);
y_opt = 0.5*[y_hl,y_out];

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
subplot(2,1,1)
hold on
plot(x,y)
scatter(VaR, 0,'MarkerFaceColor','red')
hold off
title('Profit/loss density and $$99\%$$ VaR')
set(gca,'TickLabelInterpreter','latex')

subplot(2,1,2)
hold on
plot(x,y_opt)
scatter(VaR, 0,'MarkerFaceColor','red')
hold off
title('Optimal IS candidate')

set(gca,'TickLabelInterpreter','latex')
name = 'figures/optimal_IS_cand.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')