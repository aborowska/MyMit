if plot_on
    figure(14)

    if strcmp(model,'sv')
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
        set(gcf,'defaulttextinterpreter','latex');
        
        subplot(3,1,1)
        autocorr(theta(:,1),50);
        [h1, ~, bounds1] = autocorr(theta(:,1),50);
        title('$$c$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end         
        subplot(3,1,2)
        autocorr(theta(:,2),50);
        [h2, ~, bounds2] = autocorr(theta(:,2),50);
        title('$$\phi$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(3,1,3)
        autocorr(theta(:,3),50);
        [h3, ~, bounds3] = autocorr(theta(:,3),50);
        title('$$\sigma^{2}_{\eta}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        if print_on
            name = 'figures/sv_autorr.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    else
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'defaulttextinterpreter','latex');
    
        subplot(2,2,1)
        autocorr(theta(:,1),50);
        [h1, ~, bounds1] = autocorr(theta(:,1),50);
        title('$$c$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        subplot(2,2,2)
        autocorr(theta(:,2),50);
        [h2, ~, bounds2] = autocorr(theta(:,2),50);
        title('$$\phi$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        subplot(2,2,3)
        autocorr(theta(:,3),50);
        [h3, ~, bounds3] = autocorr(theta(:,3),50);
        title('$$\sigma^{2}_{\eta}$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
        
        
        subplot(2,2,4)
        autocorr(theta(:,4),50);
        [h4, ~, bounds4] = autocorr(theta(:,4),50);
        title('$$\nu$$')
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end
       
        if print_on
            name = 'figures/svt_autorr.png';
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
    end
    L1 = min(find((h1<bounds1(1,1)) & (h1>bounds1(2,1))));
    L2 = min(find((h2<bounds2(1,1)) & (h2>bounds2(2,1))));
    L3 = min(find((h3<bounds3(1,1)) & (h3>bounds3(2,1))));
    IF1 = 1 + 2*sum(h1(1:L1));
    IF2 = 1 + 2*sum(h2(1:L2));
    IF3 = 1 + 2*sum(h3(1:L3));
    if strcmp(model,'svt')
        L4 = min(find((h4<bounds4(1,1)) & (h4>bounds4(2,1))));
        IF4 = 1 + 2*sum(h4(1:L4));
    end
end