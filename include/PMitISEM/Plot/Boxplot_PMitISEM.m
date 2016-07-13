function [VaR_outlier, ES_outlier] = Boxplot_PMitISEM(VaR_prelim, VaR_IS, ES_prelim, ES_IS, model, algo, H, N_sim, save_on, labels_in)
    close all
    %%
    ff = figure(6); 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    if (nargin==10)
        L1 = ['VaR ',labels_in{1}];
        L2 = ['VaR ',labels_in{2}];
        VB = boxplot([VaR_prelim,VaR_IS],'labels',{L1, L2});
    else     
        VB = boxplot([VaR_prelim,VaR_IS],'labels',{'VaR adapt','VaR pmit'});
    end
    set(findobj(gca,'Type','text'),'FontSize',12);
    set(VB,'LineWidth',1);
    VaR_outlier = findobj(VB,'Tag','Outliers');
    VaR_outlier = get(VaR_outlier,'YData');
    
%     GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition'); 
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.5*GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-2.5*(GT(1)-GT(4))];
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D('FontSize',12);
    clear GP GO GT move
    
    if save_on
%         name = ['figures/PMitISEM/',model,'_',algo,'_VaR_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        name = ['figures/PMitISEM/',model,'_',algo,'_VaR_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];        
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try 
%                 print(name,'-dpng','-r0')
%                 print(name,'-depsc','-r0')                    
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
        end
    end
    
    %%  
    ff = figure(7); 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    
    if (nargin==10)
        L1 = ['ES ',labels_in{1}];
        L2 = ['ES ',labels_in{2}];
        EB = boxplot([ES_prelim,ES_IS],'labels',{L1, L2});
    else     
        EB = boxplot([ES_prelim,ES_IS],'labels',{'ES adapt','ES pmit'});
    end
    set(findobj(gca,'Type','text'),'FontSize',12);  
    set(EB,'LineWidth',1);    
    ES_outlier = findobj(EB,'Tag','Outliers');
    ES_outlier = get(ES_outlier,'YData');    
    
%     GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition');
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.5*GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-2.5*(GT(1)-GT(4))];
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D('FontSize',12);
    clear GP GO GT move
    
    if save_on
%         name = ['figures/PMitISEM/',model,'_',algo,'_ES_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        name = ['figures/PMitISEM/',model,'_',algo,'_ES_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];        
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try 
%                 print(name,'-dpng','-r0')
%                 print(name,'-depsc','-r0')                    
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
            end   
    end
end