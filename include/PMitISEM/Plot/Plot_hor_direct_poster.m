function Plot_hor_direct_poster(y_H, y_T, VaR_prelim, model, save_on)
    [M,H] = size(y_H);
    
    f_pl = @(aa) 100*(exp(aa/100) - 1); % percentage profit loss
    ret = f_pl(y_H);
    ret(imag(ret)~=0) = -Inf; 
    ret = cumsum(ret,2);
    
    y_H = [y_T*ones(M,1),ret];
    clear ret
    
    ind = isfinite(y_H(:,H+1)) & (all(abs(y_H)<(-2)*(VaR_prelim(sim,1)),2)); 
    y_H = y_H(ind,:);

    ind_red = (y_H(1:1000,H+1) <= VaR_prelim(sim,1));

    figure(10)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
%     set(gcf,'defaulttextinterpreter','latex');

    hold on
    plot(0:H,y_H(~ind_red,:)','k','LineWidth',2)
    plot(0:H,y_H(ind_red,:)','r','LineWidth',2)
    plot(0:H,VaR_prelim(sim,1)*ones(1,1+H),'m','LineWidth',3) 
    hold off
    
    xlabel('Forecast horizon','FontSize', 14) % x-axis label
    ylabel('Cumulative return', 'FontSize', 14) % y-axis label

%     GP = get(gca, 'Position');
%     GO = get(gca, 'OuterPosition');
%     GT = get(gca, 'TightInset');
%     move = GO - GT*[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
%     set(gca, 'Position',move);
    XL = get(gca,'XLabel');
    set(XL,'interpreter','latex')
    YL = get(gca,'YLabel');
    set(YL,'interpreter','latex')    
    
    [tick_sort,ind_tick] = sort([VaR_prelim(sim,1), get(gca, 'YTick')]);
    new_label = get(gca, 'YTickLabel');
    set(gca, 'YTick', sort([VaR_prelim(sim,1), get(gca, 'YTick')])); 
    new_label = ['VaR';new_label];
    new_label = new_label(ind_tick,:);
    set(gca, 'YTick', tick_sort); 
    set(gca,'YTickLabel',new_label)   
    
    plotTickLatex2D('FontSize',14);
    
%     XL = get(gca,'XLabel');
%     XLp = get(XL,'Position');
%     XLp(2) = XLp(2)+move(4);
%     set(XL,'Position',XLp)
%     YL = get(gca,'YLabel');
%     YLp = get(YL,'Position');
%     YLp(1) = YLp(1)+1-move(3);
%     set(YL,'Position',YLp)
    
    if save_on
        name = ['figures/PMitISEM/',model,'_hor_direct_H', num2str(H),'_poster.png'];
        set(gcf,'PaperPositionMode','auto');
%             print(name,'-dpng','-r0')
            print(name,'-depsc','-r0')
    end
end