  figure(10) 
%     set(gcf,'units','normalized','outerposition',[0 0 0.3 0.4]);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8]);
%     set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(0:H,y_H(~ind_red,:)','k')
    plot(0:H,y_H(ind_red,:)','r')
    plot(0:H,VaR_prelim(sim,1)*ones(1,1+H),'m','LineWidth',2) 
    hold off 
    
    xlabel('Forecast horizon') % x-axis label
    ylabel('Cumulative return') % y-axis label
    
    
%     GP = get(gca, 'Position')
%     GO = get(gca, 'OuterPosition')
%     GT = get(gca, 'TightInset') 
%     GPT = GP+GT
%     GAPP = get(gca,'ActivePositionProperty')
%     XL = get(gca,'XLabel');
%     set(XL,'units','normalized')
%     XL = get(XL,'Position');
%     XL(2)
%     YL = get(gca,'YLabel');
%     set(YL,'units','normalized')
%     YL = get(YL,'Position');
%     YL(1) 
    [tick_sort,ind_tick] = sort([VaR_prelim(sim,1), get(gca, 'YTick')]);
    % set(gca, 'YTick', sort([VaR_prelim, get(gca, 'YTick')])); 
    new_label = get(gca, 'YTickLabel');
    new_label = ['VaR';new_label];
    new_label = new_label(ind_tick,:);
    set(gca, 'YTick', tick_sort); 
    set(gca,'YTickLabel',new_label)
    
    GP = get(gca, 'Position')
    GO = get(gca, 'OuterPosition')
    GT = get(gca, 'TightInset') 
    GD(1)=GP(1)-GT(1)
    GD(2)=GP(2)-GT(2)
    GD(3)=GP(1)+GT(1)+GD(1)
    GD(4)=GP(1)-GT(1)+GD(2)
    
    SH = GP(3)+GT(1)+GT(3)
    1/SH
    SV = GP(4)+GT(2)+GT(4)
    1/SV
% GT = [0.0071    0.0118    0.0036         0]
% GP = [0.1300    0.1100    0.7750    0.8150]    
    set(gca,'OuterPosition',[-0.04 0 1.1 1.06],'units','normalized')   
     1/( 0.8150+0.13) =  1.0582
aa = 1/(GP(1)+GP(4))
     set(gca,'OuterPosition',[0 0 1 1],'units','normalized')   
     
      plotTickLatex2D;

% margin = get(haxes(1),'TightInset') * [0;0;1;0] ;
% set(haxes,'OuterPosition',[-margin 0 1+margin 1])
%     set(gca,'position',[YL(1) XL(2) 1 1],'units','normalized')   
%     set(gca,'LooseInset',get(gca,'TightInset')); 
%     set(gca,'LooseInset',[-1 -1 1 1]); 
%     set(gca, 'Box', 'off');
    