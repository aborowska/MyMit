function Boxplot_Combine(model, H, N_sim, p_bar, save_on)
    % clear all
    close all
    % 
    % addpath(genpath('include/'));

%     model =  't_garch2_noS';
    % model =  't_gas';
    % model = 'arch';
    
%     H = 250;
%     N_sim = 20;
%     p_bar = 0.01;
%     save_on = 1;

    algo = 'Direct';
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_direct','ES_direct')

    algo = 'Prelim';
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')

    algo = 'MitISEM';
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    try
        load(name,'VaR_mit','ES_mit')
    catch
        VaR_mit = NaN*ones(N_sim,1);
        ES_mit = NaN*ones(N_sim,1);
    end
    algo = 'PMitISEM';
    name = ['results/PMitISEM/',model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit')

    VaR_mat = [VaR_direct, VaR_prelim, VaR_mit, VaR_pmit];
    ES_mat = [ES_direct, ES_prelim, ES_mit, ES_pmit];

    VaR_label = {'VaR naive','VaR adapt','VaR mit','VaR pmit'};
    ES_label = {'ES naive','ES adapt', 'ES mit','ES pmit'};

    ff = figure(1234);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
        VB = boxplot(VaR_mat,VaR_label);
        set(findobj(gca,'Type','text'),'FontSize',12);
        set(VB,'LineWidth',1); 
        GO = get(gca, 'OuterPosition'); 
        GT = get(gca, 'TightInset') ;
        if (H ~= 40)
            move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-2.0*(GT(1)-GT(4))];
        else
            move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-3.0*(GT(1)-GT(4))];
        end
        set(gca, 'Position',move);
        XL = findobj(gca, 'type', 'text');
        set(XL, 'Interpreter', 'latex');
        plotTickLatex2D('FontSize',12);
        clear GP GO GT move

        if save_on
    %         name = ['figures/PMitISEM/',model,'_VaR_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
            name = ['figures/PMitISEM/',model,'_VaR_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try 
%                     print(name,'-dpng','-r0')
%                     print(name,'-depsc','-r0')                    
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end


    ff = figure(6789);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
        EB = boxplot(ES_mat,ES_label);
        set(findobj(gca,'Type','text'),'FontSize',12);  
        set(EB,'LineWidth',1);

        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        if ((H == 10) || (H == 100))
            move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-2.0*(GT(1)-GT(4))];
        elseif (H == 20)
            move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-2.5*(GT(1)-GT(4))];
        else
            move = [1.5*GT(1), 1.5*GT(1), GO(3)-2*(GT(1)-GT(3)), GO(4)-3.0*(GT(1)-GT(4))];
        end
        set(gca, 'Position',move);
        XL = findobj(gca, 'type', 'text');
        set(XL, 'Interpreter', 'latex');
        plotTickLatex2D('FontSize',12);
        clear GP GO GT move

        if save_on
    %         name = ['figures/PMitISEM/',model,'_ES_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
            name = ['figures/PMitISEM/',model,'_ES_box_comb','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try 
%                     print(name,'-dpng','-r0')
%                     print(name,'-depsc','-r0')                    
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end

    % ERRORBARS   
    ff = figure(1111);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   

        hold on
        plot(ones(N_sim,1),VaR_mat(:,1),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(2*ones(N_sim,1),VaR_mat(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(3*ones(N_sim,1),VaR_mat(:,3),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(4*ones(N_sim,1),VaR_mat(:,4),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        errorbar(mean(VaR_mat),std(VaR_mat,1),'rx','Linewidth',2,'Markersize',10)     
        hold off

        set(gca,'XTick',1:4)
        set(gca,'XTickLabel', VaR_label)    

        GO = get(gca, 'OuterPosition'); 
        GT = get(gca, 'TightInset') ;
        move = [1.5*GT(1), 1.2*GT(1), GO(3)-2.2*(GT(1)-GT(3)), GO(4)-3.5*(GT(1)-GT(4))];
        set(gca, 'Position',move);

        plotTickLatex2D('FontSize',12);

        if save_on
    %         name = ['figures/PMitISEM/',model,'_VaR_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
            name = ['figures/PMitISEM/',model,'_VaR_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try 
%                     print(name,'-dpng','-r0')
%                     print(name,'-depsc','-r0')                    
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end


    ff = figure(2222);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   

        hold on
        plot(ones(N_sim,1),ES_mat(:,1),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(2*ones(N_sim,1),ES_mat(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(3*ones(N_sim,1),ES_mat(:,3),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        plot(4*ones(N_sim,1),ES_mat(:,4),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.9,0.9,0.9],'MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',10)
        errorbar(mean(ES_mat),std(ES_mat,1),'rx','Linewidth',2,'Markersize',10)     
        hold off

        set(gca,'XTick',1:4)
        set(gca,'XTickLabel', ES_label)

        GO = get(gca, 'OuterPosition'); 
        GT = get(gca, 'TightInset') ;
        move = [1.5*GT(1), 1.2*GT(1), GO(3)-2.2*(GT(1)-GT(3)), GO(4)-3.5*(GT(1)-GT(4))];
        set(gca, 'Position',move);

        plotTickLatex2D('F  ontSize',12);

        if save_on
    %         name = ['figures/PMitISEM/',model,'_ES_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.png'];
            name = ['figures/PMitISEM/',model,'_ES_errorbar','_H', num2str(H),'_Nsim',num2str(N_sim),'.eps'];      
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try 
%                     print(name,'-dpng','-r0')
%                     print(name,'-depsc','-r0')                    
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end
end