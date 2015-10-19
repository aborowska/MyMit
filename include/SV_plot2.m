if plot_on
    if strcmp(model,'sv')   
        load('SML_ibm_smooth.mat', 'theta_smooth', 'V_smooth');
    else
        load('SMLt_ibm_smooth.mat', 'theta_smooth', 'V_smooth');
    end
    
    if strcmp(model,'sv')
        figure(4)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(3,1,1)
        [f,xi] = ksdensity(theta2(:,1));
        hold on    
        histnorm(theta2(:,1),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$c$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(3,1,2)
        [f,xi] = ksdensity(theta2(:,2));
        hold on       
        histnorm(theta2(:,2),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\phi$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(3,1,3)
        [f,xi] = ksdensity(theta2(:,3));
        hold on       
        histnorm(theta2(:,3),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\sigma^{2}_{\eta}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

    %     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
        if print_on
            name = 'figures/sv_draws_param_hl.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(5)
        set(gcf,'units','normalized','outerposition',[0 0 1 0.33]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(1,2,1)
        [f,xi] = ksdensity(theta2(:,4));
        hold on    
        histnorm(theta2(:,4),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\eta_{T+1}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(1,2,2)
        [f,xi] = ksdensity(theta2(:,5));
        hold on       
        histnorm(theta2(:,5),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\varepsilon_{T+1}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

    %     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
        if print_on
            name = 'figures/sv_draws_err_hl.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         figure(8)
%         set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
%         set(gcf,'defaulttextinterpreter','latex');
%         hold on
%         plot(x_wmean1) 
%         plot(x_wmean2,'r')
%         plot(theta_smooth,'k','LineWidth',1)
%         hold off
%     %     title('Average signal path given the draws from $$q_{1,\zeta}(\theta|y)$$ and $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)$$.')
%         if v_new
%             set(gca,'TickLabelInterpreter','latex')
%         else
%             plotTickLatex2D;
%         end
% 
%         if print_on
%             name = 'figures/sv_mean_x_both.png';
%             fig = gcf;
%             fig.PaperPositionMode = 'auto';
%             print(name,'-dpng','-r0')
%         end
    
    else  
        figure(4)
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(2,2,1)
        [f,xi] = ksdensity(theta2(:,1));
        hold on    
        histnorm(theta2(:,1),20)   
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$c$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(2,2,2)
        [f,xi] = ksdensity(theta2(:,2));
        hold on       
        histnorm(theta2(:,2),20)   
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\phi$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(2,2,3)
        [f,xi] = ksdensity(theta2(:,3));
        hold on       
        histnorm(theta2(:,3),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\sigma^{2}_{\eta}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(2,2,4)
        [f,xi] = ksdensity(theta2(:,4));
        hold on       
        histnorm(theta2(:,4),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')    
        plot(xi,f,'r');
        hold off
        title('$$\nu$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
    %     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
        if print_on
            name = 'figures/svt_draws_param_hl.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(5)
        set(gcf,'units','normalized','outerposition',[0 0 1 0.33]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(1,2,1)
        [f,xi] = ksdensity(theta2(:,5));
        hold on    
        histnorm(theta2(:,5),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\eta_{T+1}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(1,2,2)
        [f,xi] = ksdensity(theta2(:,6));
        hold on       
        histnorm(theta2(:,6),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
        hold off
        title('$$\varepsilon_{T+1}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
    %     suptitle('Empirical approximate joint posterior density $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)y$$.' )
        if print_on
            name = 'figures/svt_draws_err_hl.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end   
 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         figure(8)
%         set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
%         set(gcf,'defaulttextinterpreter','latex');
%         hold on
%         plot(x_wmean1) 
%         plot(x_wmean2,'r')
%         plot(theta_smooth,'k','LineWidth',1)
%         hold off
%         title('Average signal path given the draws from $$q_{1,\zeta}(\theta|y)$$ and $$q_{2,\zeta}(\theta,\eta_{T+1},\varepsilon_{T+1}|y)$$.')
%         if v_new
%             set(gca,'TickLabelInterpreter','latex')
%         else
%             plotTickLatex2D;
%         end
%         
%         if print_on
%             name = 'figures/svt_mean_x_both.png';
%             fig = gcf;
%             fig.PaperPositionMode = 'auto';
%             print(name,'-dpng','-r0')
%         end
    end
end