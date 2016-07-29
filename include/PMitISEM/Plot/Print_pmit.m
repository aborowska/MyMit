function Print_pmit(model,model_tex,h,p_bar,N_sim)
    fname = ['results/PMitISEM/results_',model,'_pmit.tex'];
        
    if isempty(strfind(model,'_ML'))
        ML = false;
    else
        ML = true;
    end
    
    if ML
        param = '$\{\varepsilon_{1}\}$';
    elseif strcmp(model,'WN')
        param = '\sigma^{2}';
    elseif strcmp(model,'arch')
        param = '\alpha';
    else %'t_garch','t_gas'
        param = '\theta';
    end
 
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
end