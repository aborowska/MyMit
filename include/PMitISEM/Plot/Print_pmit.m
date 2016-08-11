function Print_pmit(model,h,p_bar,N_sim,estimation)
    
    [model_tex, ML] = fn_model_tex(model);

    if (nargin == 5)
        model_tex = [model_tex,' \\textbf{(',estimation,')}'];
        estimation = [estimation,'_'];  
        fname = ['results/PMitISEM/results_',model,'_',estimation,'pmit.tex']; 
    else
        fname = ['results/PMitISEM/results_',model,'_pmit.tex']; 
        estimation = '';
    end
    
    if strcmp(model,'WN')
        param = '\sigma^{2}';
    elseif strcmp(model,'arch')
        param = '\alpha';
    else %'t_garch','t_gas'
        param = '\theta';
    end
 
    if (nargin == 5)
        name = ['results/PMitISEM/',model,'_PMitISEM_',estimation,num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    else
        name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    end
    load(name)
    FID = fopen(fname, 'w+');
    
    fprintf(FID, '\\footnotesize{  \n');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
    
    fprintf(FID, '\\begin{longtable}{cccp{3.6cm}} \n');
    caption = ['\\caption{Partial mixture properties for $H=10$ in the ', model_tex,' model.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:pmits_',model,'} \\\\ \n'];
    fprintf(FID, label);
    fprintf(FID, ' Subset & Parameters & No. of components  & Weighted$^{*}$ $\\mu$ or $\\beta$$^{*}$ \\\\ \\hline \n');

    for ii = 1:10
        fprintf(FID, '%d & ' ,ii);
        if (ii == 1)
            if ML
                component = ['$\{\varepsilon_{1}\}$'];                
            else
                component = ['$\{(',param,',\varepsilon_{1})\}$'];
            end
        else
            component = ['$\{\varepsilon_{',num2str(ii),'}\}$'];
        end
        fprintf(FID, '%s & ' ,component);
        fprintf(FID, '%d & ' ,length( pmit(ii).p));
        g = sprintf('%6.4f \\, ', pmit(ii).p*pmit(ii).mu);
        fprintf(FID, '[%s]   \\\\ [1ex] \n',g);
    end

    fprintf(FID, '\\hline \n');
    fprintf(FID, ' \\multicolumn{4}{l}{\\footnotesize{$^{*}$Weighted with the mixture weights.}} \\\\ \n');             
    fprintf(FID, ' \\multicolumn{4}{l}{\\footnotesize{$^{**}$The mode $\\mu$ or the regression coefficients $\\beta$.}} \\\\ \n');             

    fprintf(FID, '\\end{longtable} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '} \n');    
    fclose(FID);
end