if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(2007,2012,T);
    plot(xx, y)
    set(gca,'XTickLabel',{'2007',' ','2008',' ','2009',' ','2010',' ','2011',' ','2012'})
	title('IBM log-returns $$y$$')   
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    
    if print_on
		name = 'figures/sv_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end