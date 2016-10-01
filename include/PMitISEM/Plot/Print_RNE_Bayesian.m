function Print_RNE_Bayesian(model, horizons, algos, M, N_sim, p_bar,results_path)
    model_tex = fn_model_tex(model);

    fname = [results_path,'results_',model,'_RNE.tex'];
    FID = fopen(fname, 'w+');

    fprintf(FID, '\\footnotesize{  \n');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
    fprintf(FID, '\\begin{longtable}{ccccccccccc}  \n');            

    caption = ['\\caption{RNE for the $99\\%%$ VaR and ES, in the ',...
    model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
    int2str(N_sim),'$ replications to obtain NSEs and RNEs.} \n'];

    fprintf(FID, caption);

    label = ['\\label{tab:res_RNE_',model,'} \\\\ \n'];
    fprintf(FID, label);

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

        RNE_direct = NaN;     
        RNE_prelim = NaN;
        RNE_ES_direct = NaN;
        RNE_ES_prelim = NaN;
        RNE_mit = NaN;
        RNE_pmit = NaN;

        fprintf(FID, '%s & & ',num2str(h));
        for algo = algos
            try
                name = [results_path,model,'_',char(algo),'_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
                load(name)                
            catch
            end

            try
                name = [results_path,model,'_',char(algo),'_',num2str(p_bar),'_H',num2str(h),'_RNE_Nsim',num2str(N_sim),'.mat'];               
                load(name)                    
            catch
            end                
        end

        fprintf(FID, '$[$%6.4f$]$  \\\\  \n', iqr(ES_pmit));
        fprintf(FID, '  & RNE & ');
        fprintf(FID, '%6.2f & ' , mean(RNE_direct));
        fprintf(FID, '%6.2f & ' , mean(RNE_prelim));
        fprintf(FID, '%6.4f & ' , mean(RNE_mit));
        fprintf(FID, '%6.2f &  &' , mean(RNE_pmit));

        fprintf(FID, '%6.2f & ' , mean(RNE_ES_direct));
        fprintf(FID, '%6.2f &  ', mean(RNE_ES_prelim));   
        RNE_ES_mit = mean(RNE_ES_direct)*var(VaR_direct)/var(VaR_mit);
        fprintf(FID, '%6.2f & ' , RNE_ES_mit);           
        RNE_ES_pmit = mean(RNE_ES_direct)*var(VaR_direct)/var(VaR_pmit);                
        fprintf(FID, '%6.2f   \\\\ [1ex] \n', RNE_ES_pmit);

    end 
    fprintf(FID, '\\hline \n');

    fprintf(FID, '  \\multicolumn{11}{l}{\\footnotesize{Missing value (--): it was not possible to generate the particular result with the corresponding algorithm.}} \\\\ \n');     
    fprintf(FID, '  \\multicolumn{11}{l}{\\footnotesize{IQR: interquantile range.}} \\\\ \n');             

    fprintf(FID, '\\end{longtable} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '\\normalsize \n');
    fclose(FID);

    Remove_NaN(fname);
end
