function [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_prelim, VaR_IS, ES_prelim, ES_IS, model, algo, H, N_sim, save_on, labels_in)
    %%
    figure(6) 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    if (nargin==10)
        L1 = ['VaR ',labels_in{1}];
        L2 = ['VaR ',labels_in{2}];
        VB = boxplot([VaR_prelim,VaR_IS],'labels',{L1, L2});
    else     
        VB = boxplot([VaR_prelim,VaR_IS],'labels',{'VaR prelim','VaR pmit'});
    end
    VaR_outlier = findobj(VB,'Tag','Outliers');
    VaR_outlier = get(VaR_outlier,'YData');
    
%     GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition'); 
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.5*GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-1.8*(GT(1)-GT(4))];
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D;
    clear GP GO GT move
    
    if save_on
        name = ['figures/PMitISEM/',model,'_',algo,'_VaR_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        set(gcf,'PaperPositionMode','auto');
        print(name,'-dpng','-r0')
    end
    
    %%  
    figure(7) 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    
    if (nargin==10)
        L1 = ['ES ',labels_in{1}];
        L2 = ['ES ',labels_in{2}];
        EB = boxplot([VaR_prelim,VaR_IS],'labels',{L1, L2});
    else     
        EB = boxplot([VaR_prelim,VaR_IS],'labels',{'ES prelim','ES pmit'});
    end
    ES_outlier = findobj(EB,'Tag','Outliers');
    ES_outlier = get(ES_outlier,'YData');    
    
%     GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition');
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.5*GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-1.8*(GT(1)-GT(4))];
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D;
    clear GP GO GT move
    
    if save_on
        name = ['figures/PMitISEM/',model,'_',algo,'_ES_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        set(gcf,'PaperPositionMode','auto');
        print(name,'-dpng','-r0')
    end
end