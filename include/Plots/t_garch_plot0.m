if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(1998,2008,T);
    plot(xx, data)
%     set(gca,'XTickLabel',{'1998',' ','1999',' ','2000',' '})
	title('S$$\&$$P 500 log-returns $$y$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/t_garch_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end