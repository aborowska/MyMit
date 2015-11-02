if plot_on
    figure(33)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, false);

    xx = 0:0.01:0.5;    n = length(xx);
    yy = -5:0.01:5;     m = length(yy);

    [XX,YY] = meshgrid(xx,yy);
    XX1 = reshape(XX,n*m,1);
    YY1 = reshape(YY,n*m,1);

    Val1 = kernel([XX1,YY1]);
    Val = reshape(Val1,m,n);
    P = surf(XX,YY,Val);
    set(P,'LineStyle','none')

    title('(arch) The high loss density $$p_{2}(\alpha,\varepsilon_{T+1})$$.')
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$')

    %%%%%%%%%%%%%%%%

    figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = 0:0.01:0.5;
    yy = -5:0.01:5;
    Mit2 = MitISEM_plot(mit2, 3, xx, yy, GamMat);
    title('(arch) Approximation to the high loss density $$q_{2,Mit}(\alpha,\varepsilon_{T+1})$$.')
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$','Position',[0.33,3.7,-17.6])
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    campos([-3,-30,150]);
    if print_on
        name = 'figures/q2_mit_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end