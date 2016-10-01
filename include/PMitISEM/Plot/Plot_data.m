function data_stats = Plot_data(y,time,save_on,figures_path)
% garch_old: time = [1998,2008]
% gas_crisis: time=[2005, 2016.5]
close all

ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.49 0.3]);
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.66 0.4]);

set(gcf,'defaulttextinterpreter','latex');

% subplot(1,3,1)   
    ax1 = axes('Position',[0.07 0.1 0.3 0.85],'Visible','off');
    xx = linspace(time(1,1),time(1,2),length(y));
    plot(ax1,xx,y,'k','LineWidth',1)
    set(gca,'XLim',time)

	%title('S$$\&$$P 500 log-returns $$y$$', 'FontSize', 12)
    plotTickLatex2D('FontSize',11);
    
% subplot(1,3,2)
    ax2 = axes('Position',[0.43 0.1 0.3 0.85],'Visible','off');
    axes(ax2) % sets ax1 to current axes
    hist(y,20)
    h = findobj(gca, 'Type','patch');
    set(h,'FaceColor',[0.2 0.2 0.2], 'EdgeColor','k')
    plotTickLatex2D('FontSize',11);
    
% subplot(1,3,3)
    ax3 = axes('Position',[0.76 0 1 1],'Visible','off');
    data_stats.T = length(y);
    data_stats.mean = mean(y);
    data_stats.median = median(y);
    data_stats.max = max(y);
    data_stats.min = min(y);    
    data_stats.std = std(y);
    data_stats.skewness = skewness(y);
    data_stats.kurtosis = kurtosis(y);

    axes(ax3) % sets ax1 to current axes
    descr = {'\textbf{Data descriptive statistics}';
    '';
    ['T: ',num2str(length(y))];
    ['Mean: ',sprintf('%6.4f',  mean(y))];
    ['Median: ',sprintf('%6.4f',  median(y))];
    ['Min.: ',sprintf('%6.4f',  min(y))];
    ['Max.: ',sprintf('%6.4f',  max(y))];
    ['St. Dev.: ',sprintf('%6.4f',  std(y))];
    ['Skewness: ',sprintf('%6.4f',  skewness(y))];
    ['Kurtosis: ',sprintf('%6.4f',  kurtosis(y))];
    };
    text(.0,0.5,descr,'FontSize',11)
    
    
    if save_on
%         name = [figures_path,model,'_hor_direct_H', num2str(H),'.png'];
        name = [figures_path,'data.eps'];
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try 
%                 print(name,'-dpng','-r0')
%                 print(name,'-depsc','-r0')                    
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
        end
    end
end