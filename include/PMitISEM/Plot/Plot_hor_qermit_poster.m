function Plot_hor_qermit_poster(y_H, y_pmit, sim, VaR_prelim, model, algo, save_on)
    figure(12)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
    
    ind_red1 = (y_H(1:1000,H+1) <= VaR_prelim(sim,1));
    ind_red2 = (y_pmit(1:1000,H+1) <= mean(VaR_prelim));
    hold on
    plot(0:H,y_H(~ind_red1,:)','k','LineWidth',2)
    plot(0:H,y_pmit(ind_red,:)','r','LineWidth',2)
    plot(0:H,mean(VaR_prelim)*ones(1,1+H),'m','LineWidth',3) 
    hold off
    xlabel('Forecast horizon','FontSize', 14) % x-axis label
    ylabel('Cumulative return','FontSize', 14) % y-axis label

    XL = get(gca,'XLabel');
    set(XL,'interpreter','latex')
    YL = get(gca,'YLabel');
    set(YL,'interpreter','latex')
    
    [tick_sort,ind_tick] = sort([mean(VaR_prelim), get(gca, 'YTick')]);
    % set(gca, 'YTick', sort([VaR_prelim, get(gca, 'YTick')])); 
    new_label = get(gca, 'YTickLabel');
    new_label = ['VaR';new_label];
    new_label = new_label(ind_tick,:);
    set(gca, 'YTick', tick_sort); 
    set(gca,'YTickLabel',new_label)
    plotTickLatex2D('FontSize',14);
    
    if save_on
        name = ['figures/PMitISEM/',model,'_',algo,'_hor_qermit_H', num2str(H),'_poster.png'];
        set(gcf,'PaperPositionMode','auto');
        print(name,'-dpng','-r0')
    end
end