if plot_on
    if strcmp(model,'sv') 
        figure(2)
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(3,1,1)  
        [f,xi] = ksdensity(theta1(:,1));
        hold on 
        histnorm(theta1(:,1),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
    %     f = normpdf(xi);
    %     plot(xi, f, 'k');
        hold off
        title('$$c$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(3,1,2)
        [f,xi] = ksdensity(theta1(:,2));
        hold on
        histnorm(theta1(:,2),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')    
        plot(xi,f,'r');
    %     f = exp(logpdf_beta((xi+1)/2));
    %     plot(xi, f, 'k');
        hold off
        title('$$\phi$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(3,1,3)
        [f,xi] = ksdensity(theta1(:,3));
        hold on
        histnorm(theta1(:,3),20)
        h = findobj(gca, 'Type','patch');
        set(h ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
        plot(xi,f,'r');
    %     f = exp(logpdf_invgamma(xi));
    %     plot(xi, f, 'k');
        hold off
        title('$$\sigma^{2}_{\eta}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

    %     suptitle('Empirical approximate posterior parameter density $$q_{1,\zeta}(\theta|y)$$.' )
        if print_on
            name = 'figures/sv_draws_param_init.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    else
        figure(2)
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'defaulttextinterpreter','latex');

        subplot(2,2,1)  
        [f,xi] = ksdensity(theta1(:,1));
        hold on
        histnorm(theta1(:,1),20)
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
        [f,xi] = ksdensity(theta1(:,2));
        hold on
        histnorm(theta1(:,2),20)
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
        [f,xi] = ksdensity(theta1(:,3));
        hold on
        histnorm(theta1(:,3),20)
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
        [f,xi] = ksdensity(theta1(:,4));
        hold on
        histnorm(theta1(:,4),20)
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
        
    %     suptitle('Empirical approximate posterior parameter density $$q_{1,\zeta}(\theta|y)y$$.' )
        if print_on
            name = 'figures/svt_draws_param_init.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
end