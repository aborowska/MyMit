if plot_on
    c = 100*log(1 + VaR_prelim/100);
    eps_bord = c./stdev;
    high_loss_subspace = [alpha1, eps_bord];
    
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = min(alpha1(:,1)):0.01:max(alpha1(:,1));
    yy = c./f_stdev(xx);
    hold on
    xlabel('$$\alpha$$')
    ylabel('$$\varepsilon_{T+1}$$')
    % scatter(high_loss_subspace(:,1),high_loss_subspace(:,2),10,'r.')
    plot(xx, yy, 'r')
    % h = area(xx,yy);
    % h.FaceColor = [1,0.4,0.3];
    % h.EdgeColor = 'red';
    scatter(alpha1(:,1),eps1(:,1),10,'k.')
    hold off
    title('(arch) Joint density $$p(\alpha,\varepsilon_{T+1}|y)$$ and hight loss subspace')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end

    if print_on
        name = 'figures/high_loss_mitisem.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end