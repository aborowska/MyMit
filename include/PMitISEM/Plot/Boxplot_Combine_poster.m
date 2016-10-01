close all
addpath(genpath('include/'));

model =  't_garch2_noS';
H = 10;
N_sim = 20;
p_bar = 0.01;
save_on = 1;

algo = 'Direct';
name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name,'VaR_direct','ES_direct')

algo = 'Prelim';
name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name,'VaR_prelim','ES_prelim')

algo = 'MitISEM';
name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name,'VaR_mit','ES_mit')

algo = 'PMitISEM';
name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name,'VaR_pmit','ES_pmit')

if H >= 100
    VaR_mit = NaN*VaR_pmit;
    ES_mit = NaN*ES_pmit;
end

VaR_mat = [VaR_direct, VaR_prelim, VaR_mit, VaR_pmit];
ES_mat = [ES_direct, ES_prelim, ES_mit, ES_pmit];

VaR_label = {'VaR naive','VaR adapt','VaR mit','VaR pmit'};
ES_label = {'ES naive','ES adapt', 'ES mit','ES pmit'};
    
figure(1234)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
    VB = boxplot(VaR_mat,VaR_label);
    set(findobj(gca,'Type','text'),'FontSize',14);
    set(VB,'LineWidth',2);
%     GO = get(gca, 'OuterPosition'); 
%     GT = get(gca, 'TightInset') ;
%     move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-2.0*(GT(1)-GT(4))];
%     set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D('FontSize',14);
%     clear GP GO GT move
    
    if save_on
        name = ['figures/PMitISEM/',model,'_VaR_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'_poster.png'];
        set(gcf,'PaperPositionMode','auto');
%             print(name,'-dpng','-r0')
            print(name,'-depsc','-r0')
    end
    
    
figure(6789)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
    EB = boxplot(ES_mat,ES_label);
    set(findobj(gca,'Type','text'),'FontSize',14);  
    set(EB,'LineWidth',2);
%     GO = get(gca, 'OuterPosition');
%     GT = get(gca, 'TightInset') ;
%     move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-2.50*(GT(1)-GT(4))];
    set(gca, 'Position',move);
    XL = findobj(gca, 'type', 'text');
    set(XL, 'Interpreter', 'latex');
    plotTickLatex2D('FontSize',14);
    clear GP GO GT move
    
    if save_on
        name = ['figures/PMitISEM/',model,'_ES_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'_poster.png'];
        set(gcf,'PaperPositionMode','auto');
%             print(name,'-dpng','-r0')
            print(name,'-depsc','-r0')
    end
    
% ERRORBARS   
figure(1111)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   

    hold on
    plot(ones(N_sim,1),VaR_mat(:,1),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(2*ones(N_sim,1),VaR_mat(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(3*ones(N_sim,1),VaR_mat(:,3),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(4*ones(N_sim,1),VaR_mat(:,4),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    errorbar(mean(VaR_mat),std(VaR_mat,1),'rx','Linewidth',2,'Markersize',10)     
    hold off
    
    set(gca,'XTick',1:4)
    set(gca,'XTickLabel', VaR_label)    
    
    GO = get(gca, 'OuterPosition'); 
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.2*GT(1), GO(3)-2.2*(GT(1)-GT(3)), GO(4)-3.5*(GT(1)-GT(4))];
    set(gca, 'Position',move);

    plotTickLatex2D;   
    
    if save_on
        name = ['figures/PMitISEM/',model,'_VaR_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        set(gcf,'PaperPositionMode','auto');
%             print(name,'-dpng','-r0')
            print(name,'-depsc','-r0')
    end
    
    
figure(2222)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   

    hold on
    plot(ones(N_sim,1),ES_mat(:,1),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(2*ones(N_sim,1),ES_mat(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(3*ones(N_sim,1),ES_mat(:,3),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    plot(4*ones(N_sim,1),ES_mat(:,4),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7])
    errorbar(mean(ES_mat),std(ES_mat,1),'rx','Linewidth',2,'Markersize',10)     
    hold off

    set(gca,'XTick',1:4)
    set(gca,'XTickLabel', ES_label)
  
    GO = get(gca, 'OuterPosition'); 
    GT = get(gca, 'TightInset') ;
    move = [1.5*GT(1), 1.2*GT(1), GO(3)-2.2*(GT(1)-GT(3)), GO(4)-3.5*(GT(1)-GT(4))];
    set(gca, 'Position',move);
     
    plotTickLatex2D;   
    
    if save_on
        name = ['figures/PMitISEM/',model,'_ES_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
        set(gcf,'PaperPositionMode','auto');
%             print(name,'-dpng','-r0')
            print(name,'-depsc','-r0')
    end