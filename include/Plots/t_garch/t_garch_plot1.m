if plot_on
    [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel, false);
    draw1 = [draw1, trnd(draw1(:,4))]; % ERRORS ARE T!!

    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'defaulttextinterpreter','latex');

    subplot(2,2,1)
    [f,xi] = ksdensity(draw1(:,1));
    plot(xi,f);
    title('$$\alpha$$')
    plotTickLatex2D;

    subplot(2,2,2)
    [f,xi] = ksdensity(draw1(:,2));
    plot(xi,f);
    title('$$\beta$$')
    plotTickLatex2D;

    subplot(2,2,3)
    [f,xi] = ksdensity(draw1(:,3));
    plot(xi,f);
    title('$$\mu$$')
    plotTickLatex2D;

    subplot(2,2,4)
    [f,xi] = ksdensity(draw1(:,4));
    plot(xi,f);
    title('$$\nu$$')
    plotTickLatex2D;

    suptitle('Smoothed kernel estimates')

    if print_on
        name = 'figures/t_garch_mitisem_draw1.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end