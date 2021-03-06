function Plot_time_precision(model, save_on, H, p_bar, N_sim, M)
    close all
    if (nargin == 2)
        H = 10;
        p_bar = 0.01;
        N_sim = 20;
        M = 10000;
    end
    
    % load results
    name = ['results/PMitISEM/',model,'_Direct_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_direct','ES_direct','time_direct')
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim','time_prelim','time_bigdraw')
    try
        name = ['results/PMitISEM/',model,'_MitISEM_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'VaR_mit','ES_mit','time_mit') 
        mit_ok = true;
    catch
        VaR_mit = NaN*ones(N_sim,1);
        ES_mit = NaN*ones(N_sim,1);
        time_mit = NaN*ones(2,1);
        mit_ok = false;
    end
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit','time_pmit')
    
    %% compute precisions
    precision_var_direct = 1/var(VaR_direct);
    precision_var_prelim = 1/var(VaR_prelim);
    precision_var_mit = 1/var(VaR_mit);
    precision_var_pmit_3 = 1/var(VaR_pmit_3);
    precision_var_pmit_4 = 1/var(VaR_pmit_4);
    precision_var_pmit_5 = 1/var(VaR_pmit_5);
    precision_var_pmit_6 = 1/var(VaR_pmit_6);
         
    precision_es_direct = 1/var(ES_direct);
    precision_es_prelim = 1/var(ES_prelim);
    precision_es_mit = 1/var(ES_mit);
    precision_es_pmit_3 = 1/var(ES_pmit_3);
    precision_es_pmit_4 = 1/var(ES_pmit_4);
    precision_es_pmit_5 = 1/var(ES_pmit_5);
    precision_es_pmit_6 = 1/var(ES_pmit_6);

    
    
    
    precision_one_digit = (1.96/0.05)^2; % 1.96*NSE<=0.05 
    
    % construct lines of time-precision tradeoff
    % construction time: time_xxx(1,1)
    % time of sampling of 10000 draws: time_xxx(2,1)
    time_direct_total = time_direct(1,1) + time_direct(2,1);
    time_prelim_total = time_prelim(1,1) + time_prelim(2,1);
    time_mit_total =  time_prelim(1,1) + time_mit(1,1) + time_mit(2,1);
    time_pmit_total_3 =  time_prelim(1,1) + time_pmit_3(1,1) + time_pmit_3(2,1);
    time_pmit_total_4 =  time_prelim(1,1) + time_pmit_4(1,1) + time_pmit_4(2,1);
    time_pmit_total_5 =  time_prelim(1,1) + time_pmit_5(1,1) + time_pmit_5(2,1);
    time_pmit_total_6 =  time_prelim(1,1) + time_pmit_6(1,1) + time_pmit_6(2,1);

    xx = 0:0.01:1500;

    slope_var_direct = precision_var_direct/time_direct(2,1);
    slope_var_prelim = precision_var_prelim/time_prelim(2,1);
    slope_var_mit = precision_var_mit/time_mit(2,1);
    slope_var_pmit_3 = precision_var_pmit_3/time_pmit_3(2,1);
    slope_var_pmit_4 = precision_var_pmit_4/time_pmit_4(2,1);
    slope_var_pmit_5 = precision_var_pmit_5/time_pmit_5(2,1);
    slope_var_pmit_6 = precision_var_pmit_6/time_pmit_6(2,1);

    slope_es_direct = precision_es_direct/time_direct(2,1);
    slope_es_prelim = precision_es_prelim/time_prelim(2,1);
    slope_es_mit = precision_es_mit/time_mit(2,1);
    slope_es_pmit_3 = precision_es_pmit_3/time_pmit_3(2,1);
    slope_es_pmit_4 = precision_es_pmit_4/time_pmit_4(2,1);
    slope_es_pmit_5 = precision_es_pmit_5/time_pmit_5(2,1);
    slope_es_pmit_6 = precision_es_pmit_6/time_pmit_6(2,1);


    line_var_direct = max(slope_var_direct*(xx - time_direct(1,1)),0);
    line_var_prelim = max(slope_var_prelim*(xx - time_prelim(1,1)),0);
    if mit_ok
        line_var_mit = max(slope_var_mit*(xx - time_prelim(1,1) - time_mit(1,1)),0);
    end
    line_var_pmit_3 = max(slope_var_pmit_3*(xx - time_prelim(1,1) - time_pmit_3(1,1)),0);
    line_var_pmit_4 = max(slope_var_pmit_4*(xx - time_prelim(1,1) - time_pmit_4(1,1)),0);
    line_var_pmit_5 = max(slope_var_pmit_5*(xx - time_prelim(1,1) - time_pmit_5(1,1)),0);
    line_var_pmit_6 = max(slope_var_pmit_6*(xx - time_prelim(1,1) - time_pmit_6(1,1)),0);

    line_es_direct = max(slope_es_direct*(xx - time_direct(1,1)),0);
    line_es_prelim = max(slope_es_prelim*(xx - time_prelim(1,1)),0);
    if mit_ok
        line_es_mit = max(slope_es_mit*(xx  - time_prelim(1,1) - time_mit(1,1)),0);
    end
    line_es_pmit_3 = max(slope_es_pmit_3*(xx - time_prelim(1,1) - time_pmit_3(1,1)),0);
    line_es_pmit_4 = max(slope_es_pmit_4*(xx - time_prelim(1,1) - time_pmit_4(1,1)),0);
    line_es_pmit_5 = max(slope_es_pmit_5*(xx - time_prelim(1,1) - time_pmit_5(1,1)),0);
    line_es_pmit_6 = max(slope_es_pmit_6*(xx - time_prelim(1,1) - time_pmit_6(1,1)),0);

    line_one_digit = precision_one_digit + 0*xx;

    %% Required for 1 digit precision 
    % time
    required_var_direct = sum(line_var_direct <= precision_one_digit)/100;
    required_var_prelim = sum(line_var_prelim <= precision_one_digit)/100;
    if mit_ok
        required_var_mit = sum(line_var_mit <= precision_one_digit)/100;
    else
        required_var_mit = NaN;
    end
    required_var_pmit_3 = sum(line_var_pmit_3 <= precision_one_digit)/100;
	required_var_pmit_4 = sum(line_var_pmit_4 <= precision_one_digit)/100;
    required_var_pmit_5 = sum(line_var_pmit_5 <= precision_one_digit)/100;
    required_var_pmit_6 = sum(line_var_pmit_6 <= precision_one_digit)/100;

    
    required_es_direct = sum(line_es_direct <= precision_one_digit)/100;
    required_es_prelim = sum(line_es_prelim <= precision_one_digit)/100;
    if mit_ok    
        required_es_mit = sum(line_es_mit <= precision_one_digit)/100;
    else
        required_es_mit = NaN;
    end
    required_es_pmit_3 = sum(line_es_pmit_3 <= precision_one_digit)/100;
    required_es_pmit_4 = sum(line_es_pmit_4 <= precision_one_digit)/100;
    required_es_pmit_5 = sum(line_es_pmit_5 <= precision_one_digit)/100;
    required_es_pmit_6 = sum(line_es_pmit_6 <= precision_one_digit)/100;
    

    if mit_ok
        required = max([required_var_direct, required_var_prelim, required_var_mit,required_var_mit,...
            required_es_direct, required_es_prelim, required_es_mit, required_es_mit]);
    else
        required = max([required_var_direct, required_var_prelim,required_var_mit,...
            required_es_direct, required_es_prelim, required_es_mit]);
    end
            
    % draws
    draws_var_direct = M*(required_var_direct - time_direct(1,1))/time_direct(2,1);
    draws_var_prelim = M*(required_var_prelim - time_prelim(1,1))/time_prelim(2,1);
    if mit_ok
        draws_var_mit = M*(required_var_mit - time_prelim(1,1) - time_mit(1,1))/time_mit(2,1);
    else
        draws_var_mit = NaN;
    end
    draws_var_pmit = M*(required_var_pmit - time_prelim(1,1) - time_pmit(1,1))/time_pmit(2,1);
    
    draws_es_direct = M*(required_es_direct - time_direct(1,1))/time_direct(2,1);
    draws_es_prelim = M*(required_es_prelim - time_prelim(1,1))/time_prelim(2,1);
    if mit_ok
        draws_es_mit = M*(required_es_mit - time_prelim(1,1) - time_mit(1,1))/time_mit(2,1);
    else
        draws_es_mit = NaN;
    end
    draws_es_pmit = M*(required_es_pmit - time_prelim(1,1) - time_pmit(1,1))/time_pmit(2,1);
    
    
    % figures
    if (required < 100)
        XTmax = 100;
        Xt = 10;
    elseif (required > 500)
        XTmax = 1500;
        Xt = 150;
    else
        XTmax = 700;
        Xt = 50;
    end
    limit = [0, XTmax, 0, 2500];

    %% VaR
    ff = figure(88);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
        %     set(gcf,'defaulttextinterpreter','latex');    
        axis(limit)
        set(gca,'XTick',0:Xt:XTmax)
        set(gca,'YTick',0:500:2500)   
        hold on
        plot(xx(1:XTmax*100),line_one_digit(1:XTmax*100),'k','LineWidth',2)
        plot(xx(1:XTmax*100),line_var_direct(1:XTmax*100),'r:','LineWidth',2)
        plot(xx(1:XTmax*100),line_var_prelim(1:XTmax*100),'m:','LineWidth',2)
        if mit_ok
            plot(xx(1:XTmax*100),line_var_mit(1:XTmax*100),'c','LineWidth',2)
        end
        hold all
        plot(xx(1:XTmax*100),line_var_pmit_3(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_var_pmit_4(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_var_pmit_5(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_var_pmit_6(1:XTmax*100),'LineWidth',2)
        hold off
        xlabel('Computing time (s)','FontSize', 12) % x-axis label
        ylabel('Precision = 1/var(VaR est.)','FontSize', 12) % y-axis label
        if mit_ok
            leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit MitISEM', 'QERMit PMitISEM');            
        else
            leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit PMitISEM');
        end
        %         set(leg,'Interpreter','latex','FontSize', 10,'Location','SouthEast');
        set(leg,'Interpreter','latex','FontSize', 12,'Location','NorthEast');

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
%         YL = get(gca,'YLabel');
%         YLp = get(YL,'Position');
%         YLp(1) = YLp(1)+1-move(3);
%         set(YL,'Position',YLp)


        if save_on
%             name = ['figures/PMitISEM/',model,'_time_precision_VaR_H', num2str(H),'.png'];
            name = ['figures/PMitISEM/',model,'_time_precision_VaR_H', num2str(H),'.eps'];
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
    
    %% ES
    ff = figure(99);
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
        %     set(gcf,'defaulttextinterpreter','latex');    
        axis(limit)
        set(gca,'XTick',0:Xt:XTmax)
        set(gca,'YTick',0:500:2500)   
        hold on
        plot(xx(1:XTmax*100),line_one_digit(1:XTmax*100),'k','LineWidth',2)
        plot(xx(1:XTmax*100),line_es_direct(1:XTmax*100),'r:','LineWidth',2)
        plot(xx(1:XTmax*100),line_es_prelim(1:XTmax*100),'m:','LineWidth',2)
        if mit_ok
            plot(xx(1:XTmax*100),line_es_mit(1:XTmax*100),'c','LineWidth',2)
        end
        hold all
        plot(xx(1:XTmax*100),line_es_pmit_3(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_es_pmit_4(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_es_pmit_5(1:XTmax*100),'LineWidth',2)
        plot(xx(1:XTmax*100),line_es_pmit_6(1:XTmax*100),'LineWidth',2)
        hold off
        xlabel('Computing time (s)','FontSize', 12) % x-axis label
        ylabel('Precision = 1/var(ES est.)','FontSize', 12) % y-axis label
        if mit_ok
            leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit MitISEM', 'QERMit PMitISEM');            
        else
            leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit PMitISEM');
        end
        set(leg,'Interpreter','latex','FontSize', 12)

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
        YL = get(gca,'YLabel');
        YLp = get(YL,'Position');
        YLp(1) = YLp(1)+1-move(3);
        set(YL,'Position',YLp)


        if save_on
%             name = ['figures/PMitISEM/',model,'_time_precision_ES_H', num2str(H),'.png'];
            name = ['figures/PMitISEM/',model,'_time_precision_ES_H', num2str(H),'.eps'];
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

% PrPaper = 1/(0.24*0.24);
% PrPaper2 = 1/(0.12*0.12);
% slope = PrPaper/5.3;
% slope2 = PrPaper2/5.4;

    %% create a tex table
    switch model
    case 't_gas'
        model_tex = 'GAS(1,1)-$t$';        
    case 't_garch2_noS'
        model_tex = 'GARCH(1,1)-$t$';
    case 'arch'
        model_tex = 'ARCH(1,1)';
    case 'WN'
        model_tex = 'White Noise';
    end

    fname = ['results/PMitISEM/results_',model,'_time_precision_H',num2str(H),'.tex'];
    FID = fopen(fname, 'w+');

    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '\\centering \n');

    caption = ['\\caption{Computation time-precision trade-off for the ',num2str(H),'-days ahead  $99\\%%$ VaR and ES evaluation in ', model_tex,' model.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:time_precision_',model,'} \n'];
    fprintf(FID, label);
    fprintf(FID, '\\begin{tabular}{lcccc}  \n');
    fprintf(FID, '  & \\multicolumn{2}{c}{Direct} & \\multicolumn{2}{c}{QERMit}  \\\\ \\cline{2-5} \n');   
    fprintf(FID, '  & Naive & Adapted & MitISEM & PMitISEM  \\\\ \\hline \n');
 
    fprintf(FID, 'Total time & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', ...
        time_direct_total, time_prelim_total, time_mit_total, time_pmit_total);
    fprintf(FID, 'Time construction candidate & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', ...
        time_direct(1,1), time_prelim(1,1), time_prelim(1,1) + time_mit(1,1), time_prelim(1,1)+time_pmit(1,1));
    fprintf(FID, ' (including $q_{1}$) &   &  (%4.2f s) & (%4.2f s) & (%4.2f s) \\\\ \n',  time_prelim(1,1),  time_prelim(1,1), time_prelim(1,1));
    fprintf(FID, ' (including $q_{2}$) &   &  & (%4.2f s) & (%4.2f s) \\\\ \n',  time_mit(1,1), time_pmit(1,1));
    fprintf(FID, '$[$Initialisation for $q_{2}$$]$&   &   & $[$%4.2f s$]$ & $[$%4.2f s$]$ \\\\ \n',  time_bigdraw, time_bigdraw);
    fprintf(FID, 'Time sampling & %4.2f s & %4.2f s & %4.2f s & %4.2f s  \\\\  \n', ...
        time_direct(2,1), time_prelim(2,1), time_mit(2,1), time_pmit(2,1));
    fprintf(FID, 'Number of draws used & %d & %d & %d & %d \\\\ \n', M, M, M, M);
    fprintf(FID, 'Time per draw & %4.2f ms & %4.2f ms & %4.2f ms & %4.2f ms \\\\ \\hline \n',...
          10*time_direct(2,1)/M, 10*time_prelim(2,1)/M, 10*time_mit(2,1)/M, 10*time_pmit(2,1)/M);
    
    fprintf(FID,'\\multicolumn{5}{c}{Required for \\%% estimate with 1 digit precision} \\\\ \\hline \n');
    fprintf(FID, 'Time: &  &  &   &  \\\\ \n');
    fprintf(FID, '\\hspace{1cm} for $VaR$ & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', required_var_direct, required_var_prelim, required_var_mit, required_var_pmit);
    fprintf(FID, '\\hspace{1cm} for $ES$ & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', required_es_direct, required_es_prelim, required_es_mit, required_es_pmit);
    fprintf(FID, 'Draws: &  &  &   &  \\\\ \n');
    fprintf(FID, '\\hspace{1cm} for $VaR$ & %s & %s  & %s  & %s  \\\\ \n',...
       num2bank(draws_var_direct), num2bank(draws_var_prelim), num2bank(draws_var_mit), num2bank(draws_var_pmit));
    fprintf(FID, '\\hspace{1cm} for $ES$ & %s & %s  & %s   & %s  \\\\ \n', ...
       num2bank(draws_es_direct), num2bank(draws_es_prelim), num2bank(draws_es_mit), num2bank(draws_es_pmit));
    
    fprintf(FID, '\\hline \n');
    
%     fprintf(FID, 'Investment: &  &  &   &  \\\\ \n');
    fprintf(FID, 'Slope: &  &  &   &  \\\\ \n');
   fprintf(FID, '\\hspace{1cm} for $VaR$ & %4.2f & %4.2f  & %4.2f  & %4.2f  \\\\ \n',...
       slope_var_direct, slope_var_prelim, slope_var_mit, slope_var_pmit);
    fprintf(FID, '\\hspace{1cm} for $ES$ & %4.2f & %4.2f  & %4.2f   & %4.2f  \\\\  \\hline \n', ...
       slope_es_direct, slope_es_prelim, slope_es_mit, slope_es_pmit);
 
    fprintf(FID, '\\end{tabular} \n');
    fprintf(FID, '\\end{table} \n');
    fprintf(FID, '} \n');
    fclose(FID);

end