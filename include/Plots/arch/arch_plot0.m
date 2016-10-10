if plot_on
    figure(100)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0:0.01:0.5;
%     kernel = @(a) posterior_arch(a, data, S, true);
    yy = kernel(xx');
    yy = yy - max(yy);
    yy = exp(yy);
    plot(xx,yy);
    xlabel('$$\alpha$$')
    % title('Posterior density $$p(\alpha)$$')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = 'figures/arch/posterior.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end   

    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0:0.01:0.5;
    yy = -5:0.01:5;
    Mit1 = MitISEM_plot(mit1, 2, xx, yy, GamMat);
    xlabel('$$\alpha$$')
    title('(arch) Approximation to the posterior density $$p(\alpha)$$')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = 'figures/q1_mit_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end