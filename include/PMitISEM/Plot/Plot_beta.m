function [Beta, Sigma, nu] = Plot_beta(pmit,model,H,save_on,version)
    SS = size(pmit,2); 
    r = size(pmit(2).mu,2);
    Beta = zeros(SS,r);
    Sigma = zeros(SS-1,1);
    nu = zeros(SS-1,1);
    for ii = 2:SS
        Beta(ii-1,:) = pmit(ii).p*pmit(ii).mu;
        Sigma(ii-1,1) = pmit(ii).p*pmit(ii).Sigma;
        nu(ii-1,1) = pmit(ii).p*pmit(ii).df';
    end
    Beta(SS,:) = Beta(SS-1,:);
    
    if ((nargin == 4) || (version == 1))    
        ff = figure(9) ;
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.5]);   
        S = stairs(2:SS+1,Beta);
        set(S(1),'Color','b','LineWidth',2)
        set(S(2),'Color','r','LineWidth',2)

        xlabel('Subset') % x-axis label

        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move = [1.8*GT(1), 1.8*GT(1), GO(3)-5*(GT(1)-GT(3)), GO(4)-10*(GT(1)-GT(4))];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        set(gca,'XTick',0.5+(2:SS))
        set(gca,'XTickLabel', int2str((2:SS)'))

        plotTickLatex2D;

        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2) + 0.1*move(4);
        set(XL,'Position',XLp)

        if (r==2)
            SL = legend('constant','sum of past returns','location','SouthWest');
            set(SL,'interpreter','latex');
        end
    
    elseif (version == 2)
        ff = figure(9);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);   
        S = stairs(2:SS+1,Beta(:,2));
        set(S(1),'Color','r','LineWidth',2)

        xlabel('Subset') % x-axis label

        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move = [1.5*GT(1), 1.6*GT(1), GO(3)-3.1*(GT(1)-GT(3)), GO(4)-3.5*(GT(1)-GT(4))];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        set(gca,'XTick',0.5+(2:2:SS))
        set(gca,'XTickLabel', int2str((2:2:SS)'))

        plotTickLatex2D;

        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2) + 0.005*move(4);
        set(XL,'Position',XLp)         
    end
    
    if save_on
%         name = ['figures/PMitISEM/',model,'_beta_H', num2str(H),'.png'];
        name = ['figures/PMitISEM/',model,'_beta_H', num2str(H),'.eps'];
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try 
%                     print(name,'-dpng','-r0')
%                     print(name,'-depsc','-r0')
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
        end
    end
        
end