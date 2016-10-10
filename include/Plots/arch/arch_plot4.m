if (plot_on && (sim == N_sim))   
    figure(8)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(sort(PL_opt))
    pos =  max(find(sort(PL_opt)<=VaR_IS(sim,1)));
    scatter(pos, VaR_IS(sim,1),'MarkerFaceColor','red')
    hold off
%     subplot(3,1,1)
%     subplot(3,1,2)
%     subplot(3,1,3)
    title('Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$.')
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end

    if print_on
        name = 'figures/arch_predict.png';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
end