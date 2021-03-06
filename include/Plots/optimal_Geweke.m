addpath(genpath('include/'));
close all

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

x = -5:0.01:5;
y = tpdf(x,5);
VaR = tinv(0.01,5);
C = tcdf(VaR,5);

y_hl = y(x<=VaR);
y_out = y(x>VaR);

% C = sum(y);
% y_hl = 0.5*y_hl/sum(y_hl);
% y_out = 0.5*y_out/sum(y_out);
y_hl = 0.5*y_hl/C;
y_out = 0.5*y_out/(1-C);

y_opt = [y_hl,y_out];

figure(1)
if v_new
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
else
   set(gcf,'units','normalized','outerposition',[0 0 0.33 0.5]);
end
set(gcf,'defaulttextinterpreter','latex');
subplot(2,1,1)
hold on
% plot(x,y)
% plot(x,y,'LineWidth',3)
plot(x,y,'LineWidth',2)
% scatter(VaR, 0,'MarkerFaceColor','red')
scatter(VaR, 0,100,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold off
% title('Profit/loss density and $$99\%$$ VaR')
% title('Profit/loss density and $$99\%$$ VaR', 'FontSize', 14)
title('Profit/loss density and $$99\%$$ VaR', 'FontSize', 12)

if v_new
    set(gca,'TickLabelInterpreter','latex')
else
%     plotTickLatex2D
%     plotTickLatex2D('FontSize',14);
    plotTickLatex2D('FontSize',12);
end


subplot(2,1,2)
hold on
% plot(x,y_opt)
% plot(x,y_opt,'LineWidth',3)
plot(x,y_opt,'LineWidth',2)
% scatter(VaR, 0,'MarkerFaceColor','red')
scatter(VaR, 0,100,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold off
% title('Optimal IS candidate')
% title('Optimal IS candidate', 'FontSize', 14)
title('Optimal IS candidate', 'FontSize', 12)

if v_new
    set(gca,'TickLabelInterpreter','latex')
else
    %     plotTickLatex2D
%     plotTickLatex2D('FontSize',14);
    plotTickLatex2D('FontSize',12);
end
% name = 'figures/other/optimal_IS_cand.png';
% name = 'figures/other/optimal_IS_cand_poster.png';
% name = 'figures/other/optimal_IS_cand_poster.eps';
name = 'figures/other/optimal_IS_cand.eps';

fig = gcf;
fig.PaperPositionMode = 'auto';
% print(name,'-dpng','-r0')
print(name,'-depsc','-r0')