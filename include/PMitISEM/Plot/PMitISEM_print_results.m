clear all 
close all


M = 10000;
N_sim = 20;
p_bar = 0.01;
model = 't_gas';
% model = 't_garch2_noS';
% model = 'arch';
% model = 'WN';

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
horizons = {'10','20','40','100','250'} ;
algos = {'Direct','Prelim','MitISEM','PMitISEM'};

VaR_mat_prelim = zeros(2,5);
VaR_mat_pmit = zeros(2,5);
ES_mat_prelim = zeros(2,5);
ES_mat_pmit = zeros(2,5);

%% VaR results PMitISEM all H
fname = ['results/PMitISEM/results_PMitISEM_',model,'.tex'];
FID = fopen(fname, 'w+');
fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
    model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:res_pmit_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{ccccccc}  \n');
fprintf(FID, ' Horizon & & $VaR_{adapt}$ & $VaR_{pmit}$ & & $ES_{adapt}$ & $ES_{pmit}$ \\\\ \\hline \n');

ii = 0;
for h = horizons
    ii = ii+1;
    VaR_prelim = NaN;
    VaR_pmit = NaN;
    ES_prelim = NaN;
    ES_pmit = NaN;
    fprintf(FID, '$%s$ & & ',char(h));
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',char(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')
        name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',char(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    try
        load(name,'VaR_pmit','ES_pmit')
    catch
        
    end
    fprintf(FID, '%6.4f & ' ,mean(VaR_prelim)); VaR_mat_prelim(1,ii) = mean(VaR_prelim);
    fprintf(FID, '%6.4f & & ',mean(VaR_pmit));  VaR_mat_pmit(1,ii) = mean(VaR_pmit);
    fprintf(FID, '%6.4f & ' ,mean(ES_prelim));  ES_mat_prelim(1,ii) = mean(ES_prelim);
    fprintf(FID, '%6.4f  \\\\ \n',mean(ES_pmit));  ES_mat_pmit(1,ii) = mean(ES_pmit);

    fprintf(FID, ' & & ');
    fprintf(FID, '(%6.4f) & ' ,std(VaR_prelim)); VaR_mat_prelim(2,ii) = std(VaR_prelim);
    fprintf(FID, '(%6.4f) & & ',std(VaR_pmit)); VaR_mat_pmit(2,ii) = std(VaR_pmit);
    fprintf(FID, '(%6.4f) & ' ,std(ES_prelim)); ES_mat_prelim(2,ii) = std(ES_prelim);
    fprintf(FID, '(%6.4f)   \\\\ [1ex] \n',std(ES_pmit));  ES_mat_pmit(2,ii) = std(ES_pmit);
end 
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\end{table} \n');
fclose(FID);


%% 
ii = 0;
for h = horizons
    ii = ii+1;
    
    VaR_prelim = NaN;
    VaR_pmit = NaN;
    ES_prelim = NaN;
    ES_pmit = NaN;   
    
    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',char(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name,'VaR_prelim','ES_prelim')
    try
        name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',char(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'VaR_pmit','ES_pmit')
    catch

    end
 
    VaR_mat_prelim(1,ii) = mean(VaR_prelim);
    VaR_mat_pmit(1,ii) = mean(VaR_pmit);
    ES_mat_prelim(1,ii) = mean(ES_prelim);
    ES_mat_pmit(1,ii) = mean(ES_pmit);

    VaR_mat_prelim(2,ii) = std(VaR_prelim);
    VaR_mat_pmit(2,ii) = std(VaR_pmit);
    ES_mat_prelim(2,ii) = std(ES_prelim);
    ES_mat_pmit(2,ii) = std(ES_pmit);
end 

hor = [10,20,40,100,250];
figure(90)
plot(hor, VaR_mat_pmit(1,:),'-or')
hold on
plot(hor, VaR_mat_pmit(1,:) + VaR_mat_pmit(2,:),'--r')
plot(hor, VaR_mat_pmit(1,:) - VaR_mat_pmit(2,:),'--r')
% hold off

plot(hor, VaR_mat_prelim(1,:),'-ob')
% hold on
plot(hor, VaR_mat_prelim(1,:) + VaR_mat_prelim(2,:),'--b')
plot(hor, VaR_mat_prelim(1,:) - VaR_mat_prelim(2,:),'--b')
hold off


figure(91)
plot(hor, ES_mat_pmit(1,:),'-or')
hold on
plot(hor, ES_mat_pmit(1,:) + ES_mat_pmit(2,:),'--r')
plot(hor, ES_mat_pmit(1,:) - ES_mat_pmit(2,:),'--r')
% hold off

plot(hor, ES_mat_prelim(1,:),'-ob')
% hold on
plot(hor, ES_mat_prelim(1,:) + ES_mat_prelim(2,:),'--b')
plot(hor, ES_mat_prelim(1,:) - ES_mat_prelim(2,:),'--b')
hold off


%% VaR results diff algo 

fname = ['results/PMitISEM/results_algos_',model,'.tex'];
FID = fopen(fname, 'w+');
fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
    model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs.} \n'];
fprintf(FID, caption);

label = ['\\label{tab:res_algos_',model,'} \n'];
fprintf(FID, label);
fprintf(FID, '\\begin{tabular}{ccccccccccc}  \n');
fprintf(FID, ' H & & $VaR_{naive}$ & $VaR_{adapt}$ & $VaR_{mit}$  & $VaR_{pmit}$ &  & $ES_{naive}$ & $ES_{adapt}$ & $ES_{mit}$ & $ES_{pmit}$ \\\\ \\hline \n');

for h = horizons
    VaR_direct = NaN;
    VaR_prelim = NaN;
    VaR_mit = NaN;
    VaR_pmit = NaN;
    ES_direct = NaN;
    ES_prelim = NaN;
    ES_mit = NaN;
    ES_pmit = NaN;
    
    fprintf(FID, '%s & & ',char(h));
    for algo = algos
        try
            name = ['results/PMitISEM/',model,'_',char(algo),'_',num2str(p_bar),'_H',char(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
            load(name)
        catch

        end
    end
    
    fprintf(FID, '%6.4f & ' ,mean(VaR_direct));
    fprintf(FID, '%6.4f & ',mean(VaR_prelim));
    fprintf(FID, '%6.4f & ' ,mean(VaR_mit));
    fprintf(FID, '%6.4f & & ',mean(VaR_pmit));
    fprintf(FID, '%6.4f & ' ,mean(ES_direct));
    fprintf(FID, '%6.4f & ' ,mean(ES_prelim));
    fprintf(FID, '%6.4f & ' ,mean(ES_mit));
    fprintf(FID, '%6.4f  \\\\ \n',mean(ES_pmit));

    fprintf(FID, ' & & ');
    fprintf(FID, '(%6.4f) & ' ,std(VaR_direct));
    fprintf(FID, '(%6.4f) & ' ,std(VaR_prelim));
    fprintf(FID, '(%6.4f) & ' ,std(VaR_mit));
    fprintf(FID, '(%6.4f) & & ',std(VaR_pmit));
    fprintf(FID, '(%6.4f) & ' ,std(ES_direct));
    fprintf(FID, '(%6.4f) & ' ,std(ES_prelim));
    fprintf(FID, '(%6.4f) & ' ,std(ES_mit));
    fprintf(FID, '(%6.4f)   \\\\ [1ex] \n',std(ES_pmit));
end 
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular} \n');
fprintf(FID, '\\raggedright \n\n'); 
fprintf(FID, '\\vspace{5pt}\\footnotesize{NaN: it was not possible to generate the particular result with the corresponding algorithm.} \n');
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
    case 't_gas'
        param = '\theta';
end
h='10';
name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
load(name)
FID = fopen(fname, 'w+');

fprintf(FID, '\\begin{table}[h] \n');
fprintf(FID, '\\centering \n');

caption = ['\\caption{Partial mixture properties for $H=10$ in the ', model_tex,' model.} \n'];
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

%% Time-Precision
save_on = true;
H = 10;
Plot_time_precision(model, save_on, H, p_bar, N_sim, M)

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
