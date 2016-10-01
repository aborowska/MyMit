function Plot_time_precision2(results, model, save_on, horizons, figures_path, estimation)

    close all
    if isempty(strfind(model,'_ML'))
        ML = false;
        no_alg = 4;
        col_styl = {'r:','m:','c','b'};
        leg = {'One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit MitISEM', 'QERMit PMitISEM'};                    
    else
        ML = true;
        no_alg = 3;
        col_styl = {'r:','c','b'};
        leg = {'One Digit Precision', 'Direct Naive', 'QERMit MitISEM', 'QERMit PMitISEM'};
    end
    
    if (nargin < 6)
        estimation = '';
    else
        model_tex = ['WN(',estimation,')'];
        estimation = [estimation,'_'];
    end
    
    precision_one_digit = (1.96/0.05)^2; % 1.96*NSE<=0.05
            
    for hh = 1:length(horizons)
        ff = figure(hh);
   
        Xmax = fn_Xmax(results.time_total(hh,:));
        Xt = ceil(Xmax/10);
        limit = [0, Xmax, 0, 2500];  
        xx = 0:0.01:Xmax;
        leg_h = leg;
      
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
   
        axis(limit)
        set(gca,'XTick',0:Xt:Xmax)
        set(gca,'XTickLabel', num2str(get(gca, 'XTick')'))
        set(gca,'YTick',0:500:2500) 
         
        hold on
        plot(xx,precision_one_digit + 0*xx,'k','LineWidth',2)
        for ii = 1:no_alg
            a = results.VaR_slope(hh,ii);
            b = results.VaR_intercept(hh,ii);
            if (isnan(a) || isnan(b))
                leg_h(ii+1) = [];
            else
                plot(xx,max(0,a*xx+b),col_styl{ii},'LineWidth',2)
            end
        end
        hold off
        xlabel('Computing time (s)','FontSize', 12) % x-axis label
        ylabel('Precision = 1/var(VaR est.)','FontSize', 12) % y-axis label
        leg_h = legend(leg_h);
        set(leg_h,'Interpreter','latex','FontSize', 10,'Location','NorthEast');
        
        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move =  GO - GT* [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        YL = get(gca,'YLabel');
        set(YL,'interpreter','latex')

        plotTickLatex2D('FontSize',12);
        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2)+50*move(4);
        set(XL,'Position',XLp)
        
        if save_on
            name = [figures_path,model,'_',estimation,'time_precision_VaR_H', num2str(horizons(hh)),'.eps'];
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try                   
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end
    end
  
    for hh = 1:length(horizons)
        ff = figure(10*hh);
   
        Xmax = fn_Xmax(results.time_total(hh,:));
        Xt = ceil(Xmax/10);
        limit = [0, Xmax, 0, 2500];  
        xx = 0:0.01:Xmax;
        leg_h = leg;

        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
   
        axis(limit)
        set(gca,'XTick',0:Xt:Xmax)
        set(gca,'XTickLabel', num2str(get(gca, 'XTick')'))
        set(gca,'YTick',0:500:2500) 
        
         
        hold on
        plot(xx,precision_one_digit + 0*xx,'k','LineWidth',2)
        for ii = 1:no_alg
            a = results.ES_slope(hh,ii);
            b = results.ES_intercept(hh,ii);
            if (isnan(a) || isnan(b))
                leg_h(ii+1) = [];
            else
                plot(xx,max(0,a*xx+b),col_styl{ii},'LineWidth',2)
            end
        end
        hold off
        xlabel('Computing time (s)','FontSize', 12) % x-axis label
        ylabel('Precision = 1/var(ES est.)','FontSize', 12) % y-axis label
        leg_h = legend(leg_h);
        set(leg_h,'Interpreter','latex','FontSize', 10,'Location','NorthEast');
        
        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move =  GO - GT* [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        YL = get(gca,'YLabel');
        set(YL,'interpreter','latex')

        plotTickLatex2D('FontSize',12);
        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2)+50*move(4);
        set(XL,'Position',XLp)
        
        if save_on
            name = [figures_path,model,'_',estimation,'time_precision_ES_H', num2str(horizons(hh)),'.eps'];
            set(gcf,'PaperPositionMode','auto');
            print_fail = 1;
            while print_fail 
                try                   
                    print(ff,name,'-depsc','-r0')
                    print_fail = 0;
                catch
                    print_fail = 1;
                end
            end
        end
    end
end
     
function Xmax = fn_Xmax(x)
    Xmax = 3.0*max(x,[],2);
    ord = floor(log(Xmax)./log(10));
    next = Xmax./10^(ord);
    next = ceil(next*2)./2;
    Xmax = next.*10.^(ord);
end