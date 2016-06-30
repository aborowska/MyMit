function Plot_time_precision(model, save_on, H, p_bar, N_sim, M)

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
    load(name,'VaR_prelim','ES_prelim','time_prelim')
    name = ['results/PMitISEM/',model,'_MitISEM_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_mit','ES_mit','time_mit') 
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_pmit','ES_pmit','time_pmit')
    
    %% compute precisions
    precision_var_direct = 1/var(VaR_direct);
    precision_var_prelim = 1/var(VaR_prelim);
    precision_var_mit = 1/var(VaR_mit);
    precision_var_pmit = 1/var(VaR_pmit);

    precision_es_direct = 1/var(ES_direct);
    precision_es_prelim = 1/var(ES_prelim);
    precision_es_mit = 1/var(ES_mit);
    precision_es_pmit = 1/var(ES_pmit);

    precision_one_digit = (1.96/0.05)^2; % 1.96*NSE<=0.05 
    
    % construct lines of time-precision tradeoff
    % construction time: time_xxx(1,1)
    % time of sampling of 10000 draws: time_xxx(2,1)
    time_direct_total = time_direct(1,1) + time_direct(2,1);
    time_prelim_total = time_prelim(1,1) + time_prelim(2,1);
    time_mit_total =  time_prelim(1,1) + time_mit(1,1) + time_mit(2,1);
    time_pmit_total =  time_prelim(1,1) + time_pmit(1,1) + time_pmit(2,1);

    xx = 0:0.01:700;

    slope_var_direct = precision_var_direct/time_direct(2,1);
    slope_var_prelim = precision_var_prelim/time_prelim(2,1);
    slope_var_mit = precision_var_mit/time_mit(2,1);
    slope_var_pmit = precision_var_pmit/time_pmit(2,1);

    slope_es_direct = precision_es_direct/time_direct(2,1);
    slope_es_prelim = precision_es_prelim/time_prelim(2,1);
    slope_es_mit = precision_es_mit/time_mit(2,1);
    slope_es_pmit = precision_es_pmit/time_pmit(2,1);


    line_var_direct = max(slope_var_direct*(xx - time_direct(1,1)),0);
    line_var_prelim = max(slope_var_prelim*(xx - time_prelim(1,1)),0);
    line_var_mit = max(slope_var_mit*(xx - time_prelim(1,1) - time_mit(1,1)),0);
    line_var_pmit = max(slope_var_pmit*(xx - time_prelim(1,1) - time_pmit(1,1)),0);

    line_es_direct = max(slope_es_direct*(xx - time_direct(1,1)),0);
    line_es_prelim = max(slope_es_prelim*(xx - time_prelim(1,1)),0);
    line_es_mit = max(slope_es_mit*(xx  - time_prelim(1,1) - time_mit(1,1)),0);
    line_es_pmit = max(slope_es_pmit*(xx - time_prelim(1,1) - time_pmit(1,1)),0);

    line_one_digit = precision_one_digit + 0*xx;

    %% Required for 1 digit precision 
    % time
    required_var_direct = sum(line_var_direct <= precision_one_digit)/100;
    required_var_prelim = sum(line_var_prelim <= precision_one_digit)/100;
    required_var_mit = sum(line_var_mit <= precision_one_digit)/100;
    required_var_pmit = sum(line_var_pmit <= precision_one_digit)/100;
    
    required_es_direct = sum(line_es_direct <= precision_one_digit)/100;
    required_es_prelim = sum(line_es_prelim <= precision_one_digit)/100;
    required_es_mit = sum(line_es_mit <= precision_one_digit)/100;
    required_es_pmit = sum(line_es_pmit <= precision_one_digit)/100;
 
    % draws
    draws_var_direct = M*(required_var_direct - time_direct(1,1))/time_direct(2,1);
    draws_var_prelim = M*(required_var_prelim - time_prelim(1,1))/time_prelim(2,1);
    draws_var_mit =  M*(required_var_mit - time_prelim(1,1) - time_mit(1,1))/time_mit(2,1);
    draws_var_pmit = M*(required_var_pmit - time_prelim(1,1) - time_pmit(1,1))/time_pmit(2,1);
    
    draws_es_direct = M*(required_es_direct - time_direct(1,1))/time_direct(2,1);
    draws_es_prelim = M*(required_es_prelim - time_prelim(1,1))/time_prelim(2,1);
    draws_es_mit = M*(required_es_mit - time_prelim(1,1) - time_mit(1,1))/time_mit(2,1);
    draws_es_pmit = M*(required_es_pmit - time_prelim(1,1) - time_pmit(1,1))/time_pmit(2,1);
    
    
    % figures
    limit = [0, 700, 0, 2500];
    
    %% VaR
    figure(888) 
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
        %     set(gcf,'defaulttextinterpreter','latex');    
        axis(limit)
        set(gca,'XTick',0:50:700)
        set(gca,'YTick',0:500:2500)   
        hold on
        plot(xx,line_one_digit,'k')
        plot(xx,line_var_direct,'r:')
        plot(xx,line_var_prelim,'m:')
        plot(xx,line_var_mit,'c')
        plot(xx,line_var_pmit,'b')
        hold off
        xlabel('Computing time (s)') % x-axis label
        ylabel('Precision = 1/var(VaR est.)') % y-axis label
%         leg = legend('one digit', 'direct', 'prelim', 'mit', 'pmit');
        leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit MitISEM', 'QERMit PMitISEM');
        set(leg,'Interpreter','latex','FontSize', 10,'Location','SouthEast');

        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move =  GO - GT* [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        YL = get(gca,'YLabel');
        set(YL,'interpreter','latex')

        plotTickLatex2D;
        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2)+50*move(4);
        set(XL,'Position',XLp)
%         YL = get(gca,'YLabel');
%         YLp = get(YL,'Position');
%         YLp(1) = YLp(1)+1-move(3);
%         set(YL,'Position',YLp)


        if save_on
            name = ['figures/PMitISEM/',model,'_time_precision_VaR_H', num2str(H),'.png'];
            set(gcf,'PaperPositionMode','auto');
            print(name,'-dpng','-r0')
        end
    
    %% ES
    figure(99) 
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.3 0.4]);
        %     set(gcf,'defaulttextinterpreter','latex');    
        axis(limit)
        set(gca,'XTick',0:50:500)
        set(gca,'YTick',0:500:2500)   
        hold on
        plot(xx,line_one_digit,'k')
        plot(xx,line_es_direct,'r:')
        plot(xx,line_es_prelim,'m:')
        plot(xx,line_es_mit,'c')
        plot(xx,line_es_pmit,'b')
        hold off
        xlabel('Computing time (s)') % x-axis label
        ylabel('Precision = 1/var(ES est.)') % y-axis label
%         leg = legend('one digit', 'direct', 'prelim', 'mit', 'pmit');
        leg = legend('One Digit Precision', 'Direct Naive', 'Direct Adaptive', 'QERMit MitISEM', 'QERMit PMitISEM');
        set(leg,'Interpreter','latex','FontSize', 10)

        GO = get(gca, 'OuterPosition');
        GT = get(gca, 'TightInset') ;
        move =  GO - GT* [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
        set(gca, 'Position',move);
        XL = get(gca,'XLabel');
        set(XL,'interpreter','latex')
        YL = get(gca,'YLabel');
        set(YL,'interpreter','latex')

        plotTickLatex2D;
        XL = get(gca,'XLabel');
        XLp = get(XL,'Position');
        XLp(2) = XLp(2)+50*move(4);
        set(XL,'Position',XLp)
        YL = get(gca,'YLabel');
        YLp = get(YL,'Position');
        YLp(1) = YLp(1)+1-move(3);
        set(YL,'Position',YLp)


        if save_on
            name = ['figures/PMitISEM/',model,'_time_precision_ES_H', num2str(H),'.png'];
            set(gcf,'PaperPositionMode','auto');
            print(name,'-dpng','-r0')
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

    fname = ['results/PMitISEM/time_precision_',model,'.tex'];
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
    fprintf(FID, ' (including $q_{2}$) &   &  & (%4.2f s) & (%4.2f s) \\\\ \n',  time_mit(1,1), time_pmit(1,1));
    fprintf(FID, 'Time sampling & %4.2f s & %4.2f s & %4.2f s & %4.2f s  \\\\  \n', ...
        time_direct(2,1), time_prelim(2,1), time_mit(2,1), time_pmit(2,1));
    fprintf(FID, 'Number of draws used & %d & %d & %d & %d \\\\ \\hline \n', M, M, M, M);
%     fprintf(FID, 'Time per draw & %4.2f ms & %4.2f ms & %4.2f ms & %4.2f ms \\\\ \\hline \n',...
%          10*time_direct_total/M, 10*time_prelim_total/M, 10*time_mit_total/M, 10*time_pmit_total/M);
    
    fprintf(FID,'\\multicolumn{5}{c}{Required for \\%% estimate with 1 digit precision} \\\\ \\hline \n');
    fprintf(FID, 'Time: &  &  &   &  \\\\ \n');
    fprintf(FID, '\\hspace{1cm} for $VaR$ & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', required_var_direct, required_var_prelim, required_var_mit, required_var_pmit);
    fprintf(FID, '\\hspace{1cm} for $ES$ & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n', required_es_direct, required_es_prelim, required_es_mit, required_es_pmit);
    fprintf(FID, 'Draws: &  &  &   &  \\\\ \n');
    fprintf(FID, '\\hspace{1cm} for $VaR$ & %6.0f & %6.0f  & %6.0f  & %6.0f  \\\\ \n',...
       draws_var_direct, draws_var_prelim, draws_var_mit, draws_var_pmit);
    fprintf(FID, '\\hspace{1cm} for $ES$ & %6.0f & %6.0f  & %6.0f   & %6.0f  \\\\ \n', ...
       draws_es_direct, draws_es_prelim, draws_es_mit, draws_es_pmit);
    
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\end{tabular} \n');
    fprintf(FID, '\\end{table} \n');
    fprintf(FID, '} \n');
    fclose(FID);

end