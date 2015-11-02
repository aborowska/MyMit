if plot_on  
        if N_sim == 1
            sim = 1;
            figure(9)
        else
            figure(sim)
        end
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
        set(gcf,'defaulttextinterpreter','latex');
        hold on
        plot(sort(PL_opt),'b')
        pos =  max(find(PL_opt_h1<=VaR_IS(sim,1)));
        scatter(pos, VaR_IS(sim,1),'MarkerEdgeColor','red','MarkerFaceColor','red')
        pos =  max(find(PL_opt_h1<=VaR_prelim));
        scatter(pos, VaR_prelim,'MarkerEdgeColor','green','MarkerFaceColor','green')
        hold off

        title(['Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$. Model: ',model,'.'])
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        if print_on
            if strcmp(model,'sv')
                name = 'figures/sv_predict.png';
            else
                name = 'figures/svt_predict.png';
            end
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
end