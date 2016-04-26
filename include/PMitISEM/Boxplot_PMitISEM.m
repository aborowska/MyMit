function Boxplot_PMitISEM(VaR_prelim,VaR_IS,ES_prelim,ES_IS,model,H,N_sim)
    figure(6) 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    boxplot([VaR_prelim,VaR_IS],'labels',{'VaR prelim','VaR full pmit'})
    GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition');
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-1.5*(GT(1)-GT(4))]
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D;
    clear GP GO GT move
    
    name = ['figures/PMitISEM/',model,'_VaR_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
    set(gcf,'PaperPositionMode','auto');
    print(name,'-dpng','-r0')
    
    figure(7) 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    boxplot([ES_prelim,ES_IS],'labels',{'ES prelim','ES full pmit'})
    GP = get(gca, 'Position');
    GO = get(gca, 'OuterPosition');
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), GT(1), GO(3)-1.8*(GT(1)-GT(3)), GO(4)-1.5*(GT(1)-GT(4))]
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D;
    clear GP GO GT move

    name = ['figures/PMitISEM/',model,'_ES_box','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
    set(gcf,'PaperPositionMode','auto');
    print(name,'-dpng','-r0')
end