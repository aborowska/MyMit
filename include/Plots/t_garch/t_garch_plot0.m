if plot_on
    figure(100)

    addpath(genpath('include/'));
    y = csvread('GSPC_ret_tgarch.csv');
    y = 100*y;

    v_new = ver('symbolic');
    v_new = v_new.Release;
    if strcmp(v_new,'(R2014b)')
        v_new = 1;
    else
        v_new = 0;
    end

    if v_new
        set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
    end
    set(gcf,'defaulttextinterpreter','latex');
    
    xx = linspace(1998,2008,length(y));
%     plot(xx, y);
    plot(xx, y,'LineWidth',2)
%     set(gca,'XTickLabel',{'1998',' ','1999',' ','2000',' '})
	title('S$$\&$$P 500 log-returns $$y$$', 'FontSize', 14)
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
    %     plotTickLatex2D
        plotTickLatex2D('FontSize',14);
    end
    
	if print_on
% 		name = 'figures/t_garch_data.png';
		name = 'figures/t_garch_data_poster.png';
        fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end