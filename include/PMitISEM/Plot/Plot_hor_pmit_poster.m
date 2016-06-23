function Plot_hor_pmit_poster(y_pmit,y_T, VaR_prelim, model, algo, save_on)
    [M,H] = size(y_pmit);
    
    f_pl = @(aa) 100*(exp(aa/100) - 1); % percentage profit loss
    ret = f_pl(y_pmit);
    ret(imag(ret)~=0) = -Inf; 
    ret = cumsum(ret,2);
    
    y_pmit = [y_T*ones(M,1),ret];
    clear ret
   
    ind = isfinite(y_pmit(:,H+1)) & (all(abs(y_pmit)<(-2)*mean(VaR_prelim),2)); 
    y_pmit = y_pmit(ind,:);
    ind_red = (y_pmit(1:1000,H+1) <= mean(VaR_prelim));

    figure(11)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
%     set(gcf,'defaulttextinterpreter','latex');

    hold on
    plot(0:H,y_pmit(~ind_red,:)','k','LineWidth',2)
    plot(0:H,y_pmit(ind_red,:)','r','LineWidth',2)
    plot(0:H,mean(VaR_prelim)*ones(1,1+H),'m','LineWidth',3) 
    hold off
    xlabel('Forecast horizon','FontSize', 14) % x-axis label
    ylabel('Cumulative return','FontSize', 14) % y-axis label
    
%     GP = get(gca, 'Position');
%     GO = get(gca, 'OuterPosition');
%     GT = get(gca, 'TightInset') ;
%     move =  GO - GT* [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
%     set(gca, 'Position',move);
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

%     XL = get(gca,'XLabel');
%     XLp = get(XL,'Position');
%     XLp(2) = XLp(2)+move(4);
%     set(XL,'Position',XLp)
%     YL = get(gca,'YLabel');
%     YLp = get(YL,'Position');
%     YLp(1) = YLp(1)+1-move(3);
%     set(YL,'Position',YLp)
    
    if save_on
        name = ['figures/PMitISEM/',model,'_',algo,'_hor_pmit_H', num2str(H),'_poster.png'];
        set(gcf,'PaperPositionMode','auto');
        print(name,'-dpng','-r0')
    end
end