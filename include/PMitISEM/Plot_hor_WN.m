y_H = [y(end)*ones(M,1),y_H];

ret = cumsum(y_H,2);
ind_red = (ret(1:500,H+1) <= VaR_prelim);

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.4]);
set(gcf,'defaulttextinterpreter','latex');

hold on
plot(0:H,ret(~ind_red,:)','k')
plot(0:H,VaR_prelim*ones(1,H+1),'m','LineWidth',2) 
plot(0:H,ret(ind_red,:)','r')
hold off
xlabel('Forecast horizon') % x-axis label
ylabel('Cumulative return') % y-axis label

[tick_sort,ind_tick] = sort([VaR_prelim, get(gca, 'YTick')]);
% set(gca, 'YTick', sort([VaR_prelim, get(gca, 'YTick')])); 
new_label = get(gca, 'YTickLabel');
new_label = ['VaR';new_label];
new_label = new_label(ind_tick,:);
set(gca, 'YTick', tick_sort); 
set(gca,'YTickLabel',new_label)
plotTickLatex2D;
set(gca,'position',[0 0 1 1],'units','normalized')