function Print_PMitISEM_alg_comb(model, horizons, algos, M, N_sim, p_bar)
    [model_tex, ML] = fn_model_tex(model);
    
    if strcmp(model,'WN_ML')
        estimation = {'_true','_mle'};
    else
        estimation = {''};
    end
    
    for est = estimation       
        fname = ['results/PMitISEM/results_',model,char(est),'_alg_comp.tex'];
        FID = fopen(fname, 'w+');
        
        fprintf(FID, '\\footnotesize{  \n');
        fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
        if ML
            fprintf(FID, '\\begin{longtable}{ccccccccc}  \n');
        else
            fprintf(FID, '\\begin{longtable}{ccccccccccc}  \n');            
        end
        if strcmp(model,'WN_ML')
            est2 = char(est); 
            if strcmp(est2,'_mle')
                est2 = 'MLE';
            else
                est2 = 'true';
            end
            caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
            model_tex,' model, given the ',est2,' parameter, based on $N=',int2str(M),'$ candidate draws and $',...
            int2str(N_sim),'$ replications to obtain NSEs.} \n'];
        else
            caption = ['\\caption{Results for the $99\\%%$ VaR and ES, in the ',...
            model_tex,' model, based on $N=',int2str(M),'$ candidate draws and $',...
            int2str(N_sim),'$ replications to obtain NSEs.} \n'];
        end
        fprintf(FID, caption);

        label = ['\\label{tab:res_algos_',model,'} \\\\ \n'];
        fprintf(FID, label);
        if ML
            fprintf(FID, ' H & & $VaR_{naive}$ & $VaR_{mit}$ & $VaR_{pmit}$ &  & $ES_{naive}$ & $ES_{mit}$ & $ES_{pmit}$ \\\\ \\hline \n');
        else
            fprintf(FID, ' H & & $VaR_{naive}$ & $VaR_{adapt}$ & $VaR_{mit}$  & $VaR_{pmit}$ &  & $ES_{naive}$ & $ES_{adapt}$ & $ES_{mit}$ & $ES_{pmit}$ \\\\ \\hline \n');
        end

        for h = horizons
            VaR_direct = NaN;
            if ~ML
                VaR_prelim = NaN;
            end
            VaR_mit = NaN;
            VaR_pmit = NaN;
            ES_direct = NaN;
            if ~ML
                ES_prelim = NaN;
            end
            ES_mit = NaN;
            ES_pmit = NaN;
            
            RNE_direct = NaN;
            if ~ML
                RNE_prelim = NaN;
            end
            RNE_mit = NaN;
            RNE_pmit = NaN;
            
            fprintf(FID, '%s & & ',num2str(h));
            for algo = algos
                try
                    name = ['results/PMitISEM/',model,'_',char(algo),char(est),'_',num2str(p_bar),'_H',num2str(h),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
                    load(name)
                catch

                end
            end

            fprintf(FID, '%6.4f & ' , mean(VaR_direct));
            if ~ML
                fprintf(FID, '%6.4f & ', mean(VaR_prelim));
            end
            fprintf(FID, '%6.4f & ' , mean(VaR_mit));
            fprintf(FID, '%6.4f & & ', mean(VaR_pmit));
            fprintf(FID, '%6.4f & ' , mean(ES_direct));
            if ~ML
                fprintf(FID, '%6.4f & ' , mean(ES_prelim));
            end
            fprintf(FID, '%6.4f & ' , mean(ES_mit));
            fprintf(FID, '%6.4f  \\\\ \n', mean(ES_pmit));

            fprintf(FID, '  & NSE & ');
            fprintf(FID, '(%6.4f) & ' , std(VaR_direct));
            if ~ML
                fprintf(FID, '(%6.4f) & ' , std(VaR_prelim));      
            end
            fprintf(FID, '(%6.4f) & ' , std(VaR_mit));
            fprintf(FID, '(%6.4f) & & ', std(VaR_pmit));
            fprintf(FID, '(%6.4f) & ' , std(ES_direct));
            if ~ML
                fprintf(FID, '(%6.4f) & ' , std(ES_prelim));   
            end
            fprintf(FID, '(%6.4f) & ' , std(ES_mit));
        %     fprintf(FID, '(%6.4f)   \\\\ [1ex] \n', std(ES_pmit));
            fprintf(FID, '(%6.4f)   \\\\ \n', std(ES_pmit));

            fprintf(FID, ' & IQR & ');
            fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_direct));
            if ~ML
                 fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_prelim));       
            end
            fprintf(FID, '$[$%6.4f$]$ & ' , iqr(VaR_mit));
            fprintf(FID, '$[$%6.4f$]$ & & ', iqr(VaR_pmit));
            fprintf(FID, '$[$%6.4f$]$ & ' , iqr(ES_direct));
            if ~ML     
                fprintf(FID, '$[$%6.4f$]$  &' , iqr(ES_prelim));
            end
            fprintf(FID, '$[$%6.4f$]$ & ' , iqr(ES_mit));    
            if (~ML)
                fprintf(FID, '$[$%6.4f$]$  \\\\  \n', iqr(ES_pmit));
                fprintf(FID, '  & RNE & ');
                fprintf(FID, '%6.2f & ' , mean(RNE_direct));
                fprintf(FID, '%6.2f & ' , mean(RNE_prelim));
                fprintf(FID, '%6.2f & ' , mean(RNE_mit));
                fprintf(FID, '%6.2f &  &' , mean(RNE_pmit));

                fprintf(FID, '? %6.2f & ' , (1/std(ES_direct))^2);
                fprintf(FID, '? %6.2f &  ', (1/std(ES_prelim))^2);
                fprintf(FID, '? %6.2f & ' , (1/std(ES_mit))^2);
                fprintf(FID, '? %6.2f   \\\\ [1ex] \n', (1/std(ES_pmit))^2);
            else
                fprintf(FID, '$[$%6.4f$]$  \\\\  \n', iqr(ES_pmit));
                fprintf(FID, '  & RNE & ');
                fprintf(FID, '%6.2f & ' , 1);
                fprintf(FID, '%6.2f & ' , (1/std(VaR_mit))^2);
                fprintf(FID, '%6.2f & & ', (1/std(VaR_pmit))^2);
                fprintf(FID, '%6.2f & ' , 1);
                fprintf(FID, '%6.2f & ' , (1/std(ES_mit))^2);
                fprintf(FID, '%6.2f   \\\\ [1ex] \n', (1/std(ES_pmit))^2);
            end

        end 
        fprintf(FID, '\\hline \n');
        if ML
            fprintf(FID, '  \\multicolumn{9}{l}{\\footnotesize{NaN: it was not possible to generate the particular result with the corresponding algorithm.}} \\\\ \n');     
            fprintf(FID, '  \\multicolumn{9}{l}{\\footnotesize{IQR: interquantile range.}} \\\\ \n');                 
        else
            fprintf(FID, '  \\multicolumn{11}{l}{\\footnotesize{NaN: it was not possible to generate the particular result with the corresponding algorithm.}} \\\\ \n');     
            fprintf(FID, '  \\multicolumn{11}{l}{\\footnotesize{IQR: interquantile range.}} \\\\ \n');             
        end
        fprintf(FID, '\\end{longtable} \n');
        fprintf(FID, '} \n');
        fprintf(FID, '} \n');
        fclose(FID);
    end
end