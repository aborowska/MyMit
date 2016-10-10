if plot_on  
        if N_sim == 1
            sim = 1;
            figure(9)
        else
            figure(sim)
        end
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'defaulttextinterpreter','latex');

f_pl = @(aa) 100*(exp(aa/100) - 1); 

PL_opt_h1 = f_pl(sum(y_opt_h1,2));
PL_opt_h1(imag(PL_opt_h1)~=0) = -Inf;
[PL_opt_h1, ind] = sort(PL_opt_h1); 

hold on
plot(PL_opt_h1,'b')
pos =  max(find(PL_opt_h1<=VaR_IS(sim,1)));
scatter(pos, mean(VaR_IS),'MarkerEdgeColor','red','MarkerFaceColor','red')
pos =  max(find(PL_opt_h1<=VaR_prelim));
scatter(pos, VaR_prelim,'MarkerEdgeColor','green','MarkerFaceColor','green')
hold off

        title(['Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$. Model: $$',model,'$$.'])
        if v_new
            set(gca,'TickLabelInterpreter','latex')
        else
            plotTickLatex2D;
        end

        if print_on
            if strcmp(model,'sv_x')
                name = 'figures/sv_predict.png';
            else
                name = 'figures/svt_predict.png';
            end
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(name,'-dpng','-r0')
        end
end