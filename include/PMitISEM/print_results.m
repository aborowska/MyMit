fname = ['results/PMitISEM/results_',model,'.tex'];
FID = fopen(fname, 'w+');
fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Results for the 1-day-ahead $99\\%%$ VaR and ES, in the ',...
    model,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:res_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{ccccccc}  \n');
fprintf(FID, ' Horizon & & $VaR_{prelim}$ & $VaR_{IS}$ & & $ES_{prelim}$ & $ES_{IS}$ \\\\ \\hline \n');

horizons = {'10','100','250'} ;
for  ii=1:3
    h = char(horizons(ii));
    load(['results/PMitISEM/',model,'_0.01_H',h,'_VaR_results_Nsim20.mat'])
    fprintf(FID, '$%s$ & & ',h);
    fprintf(FID, '%6.4f & ' ,mean(VaR_prelim));
    fprintf(FID, '%6.4f & & ',mean(VaR_IS));
    fprintf(FID, '%6.4f & ' ,mean(ES_prelim));
    fprintf(FID, '%6.4f  \\\\ \n',mean(ES_IS));

    fprintf(FID, ' & & ');
    fprintf(FID, '(%6.4f) & ' ,std(VaR_prelim));
    fprintf(FID, '(%6.4f) & & ',std(VaR_IS));
    fprintf(FID, '(%6.4f) & ' ,std(ES_prelim));
    fprintf(FID, '(%6.4f)   \\\\ [1ex] \n',std(ES_IS));
end 
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);


%% Box Plots
for ii=1:3
    h = char(horizons(ii));
    load(['results/PMitISEM/',model,'_0.01_H',h,'_VaR_results_Nsim20.mat'])

    figure(ii*6) 
    set(gcf,'units','normalized','outerposition',[0 0 0.3 0.4]);   
    set(gcf,'defaulttextinterpreter','latex');
    boxplot([VaR_prelim,VaR_IS],'labels',{'VaR prelim','VaR PMitISEM'})
    lab = findobj(gca, 'type', 'text');
    set(lab, 'Interpreter', 'latex');
    plotTickLatex2D;
    set(gca,'position',[0 0 1 1],'units','normalized')
    name = ['figures/PMitISEM/',model,'_',num2str(p_bar),'_H', h,'_VaR_box_Nsim',num2str(N_sim),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
 
    
    figure(ii*7) 
    set(gcf,'units','normalized','outerposition',[0 0 0.3 0.4]);   
    set(gcf,'defaulttextinterpreter','latex');
    boxplot([ES_prelim,ES_IS],'labels',{'ES prelim','ES PMitISEM'})
    lab = findobj(gca, 'type', 'text');
    set(lab, 'Interpreter', 'latex');
    plotTickLatex2D;
    set(gca,'position',[0 0 1 1],'units','normalized')
    name = ['figures/PMitISEM/',model,'_',num2str(p_bar),'_H', h,'_ES_box_Nsim',num2str(N_sim),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end
