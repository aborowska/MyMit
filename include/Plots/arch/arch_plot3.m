if (plot_on && (sim == N_sim))
    figure(4)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    scatter(draw_opt(:,1),draw_opt(:,2),10,'k.')
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$')
    hold off
    title('(arch) Draws from the optimal importance density $$q_{opt}(\alpha,\varepsilon_{T+1}|y)$$.')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = 'figures/q_opt_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(5)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0:0.01:0.5; xx = xx';
    yy = -5:0.01:5;
    n = length(xx);
    mit_ep = normpdf(yy);
    Mit1 = diag(Mit1)*repmat(mit_ep,n,1); Mit1 = Mit1';
    [X1,X2] = meshgrid(xx,yy);
    MP = surf(X1,X2,Mit1);
    set(MP,'LineStyle','none')
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$','Position',[0.57,6.3,-3.6])
    title('(arch) Joint density $$p(\alpha,\varepsilon_{T+1}|y)$$.')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    campos([-3,-32,20]);

    if print_on
        name = 'figures/joint_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 

    figure(6)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    Mit_opt = 0.5*Mit1 + 0.5*Mit2;
    surf(X1,X2,Mit_opt,'EdgeColor','interp'); %colormap(bone); colormap(hot)
    title('(arch) Approximation to the optimal posterior density')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    campos([-3,-32,60]);
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$','Position',[0.5,5.5,-10])

    if print_on
        name = 'figures/q_opt_mit_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end    
end