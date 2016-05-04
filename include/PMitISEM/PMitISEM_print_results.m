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
    VaR_prelim = NaN;
    VaR_IS = NaN;
    ES_prelim = NaN;
    ES_IS = NaN;
    h = char(horizons(ii));
    fprintf(FID, '$%s$ & & ',h);
    try
        load(['results/PMitISEM/',model,'_0.01_H',h,'_VaR_results_Nsim20.mat'],'VaR_prelim','VaR_IS','ES_prelim','ES_IS')
    catch

    end
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



%% Pmit print
fname = ['results/PMitISEM/pmit_',model,'.tex'];
switch model
    case 'WN'
        param = '\sigma^{2}';
    case 'arch'
        param = '\alpha';
    case 't_garch'
        param = '\theta';
end
h='10';
load(['results/PMitISEM/',model,'_0.01_H',h,'_VaR_results_Nsim20.mat'],'pmit')
FID = fopen(fname, 'w+');

fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Partial mixture properties for $H=10$ in the ', model,' model.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:pmits_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{cccc}  \n');
fprintf(FID, ' Subset & Parameters& No. of components $h_{s}$ & weighted $\\mu$ or $\\beta$  \\\\ \\hline \n');

for ii = 1:10
    fprintf(FID, '%d & ' ,ii);
    if (ii == 1)
        component = ['$\{(',param,',\varepsilon_{1})\}$'];
    else
        component = ['$\{\varepsilon_{',num2str(ii),'}\}$'];
    end
	fprintf(FID, '%s & ' ,component);
    fprintf(FID, '%d & ' ,length( pmit(ii).p));
    g = sprintf('%6.4f, ', pmit(ii).p*pmit(ii).mu);
    fprintf(FID, '[%s]   \\\\ [1ex] \n',g);
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
