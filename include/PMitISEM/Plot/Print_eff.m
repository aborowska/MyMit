function Print_eff(model, horizons,results_path)
  
    M = 10000;
    BurnIn = 1000;
    N_sim = 20;
    p_bar = 0.01;
 
    [model_tex, ML] = fn_model_tex(model);
 
    %% Latex table
    fname = [results_path,'results_',model,'_eff.tex'];

    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{longtable}{cc cc c cc} \n');

    caption = ['\\caption{Efficiency of the high-loss space approximations in ', model_tex,...
        ' model:\\\\ basic MitISEM vs. PMitISEM.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:eff',model,'} \\\\ \n'];
    fprintf(FID, label);
    fprintf(FID, ' && \\multicolumn{2}{c}{MitISEM} && \\multicolumn{2}{c}{PMitISEM} \\\\ \\cline{3-4} \\cline{6-7} \n'); 
    fprintf(FID, ' $H$ && Eff  & RNE && Eff & RNE  \\\\ \\cline{1-1} \\cline{3-4} \\cline{6-7} \n'); 

    for hh = 1: length(horizons)
        fprintf(FID, '%s & &', num2str(horizons(hh)));
        
        algo = 'MitISEM';
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(horizons(hh)),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        try
            load(name,'mit_eff')
            fprintf(FID, '%6.4f & ', mit_eff );
        catch
            fprintf(FID, '-- & ');       
        end
        try
            load(name,'RNE_mit')
            fprintf(FID, ' %6.4f && ', mean(RNE_mit) );
        catch
            fprintf(FID, '-- && ');       
        end
        algo = 'PMitISEM';
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(horizons(hh)),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'RNE_pmit','pmit_eff')
        fprintf(FID, '%6.4f &  %6.4f \\\\ [1ex] \n', pmit_eff, mean(RNE_pmit) );  
    end
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\multicolumn{7}{p{7cm}}{\\footnotesize{Eff: the fraction of draws generated from an approximation that result in high-losses.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{7}{l}{\\footnotesize{RNE: for Bayesian $VaR$ estimation.}} \n');
  
    fprintf(FID, '\\end{longtable} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '\\normalsize \n');
    
    fclose(FID);
end