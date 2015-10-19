if strcmp(s,'direct')
    fileID = fopen(['results/',model,'_',s,'_results_',num2str(hp),'_df_',num2str(df),'.txt'],'w');
else
    fileID = fopen(['results/',model,'_',s,'_results_',num2str(hp),'.txt'],'w');
end

if strcmp(s,'direct')
    fprintf(fileID,'Controls: \n');
    fprintf(fileID,'----------------------------------------\n');
    fprintf(fileID,'df:  %2i\n',df);
    fprintf(fileID,'\n');
    fprintf(fileID,'Direct: \n');
    fprintf(fileID,'----------------------------------------\n');
    fprintf(fileID,'Direct MH acceptance rate: %6.4f\n', accept_direct);
    fprintf(fileID,'Direct (mean) VAR estimate: %6.4f\n',mean_VaR_direct);
    fprintf(fileID,'Direct (mean) ES estimate: %6.4f\n',mean_ES_direct);
    fprintf(fileID,'\n');
    fprintf(fileID,'NSE: \n');
    fprintf(fileID,'----------------------------------------\n');
    if exist('NSE_VaR_direct','var')
        fprintf(fileID,'NSE direct VaR estimate: %6.4f.\n',NSE_VaR_direct);
    end
    if exist('NSE_ES_direct','var')
        fprintf(fileID,'NSE direct ES estimate: %6.4f.\n',NSE_ES_direct);
    end
else % mit
    fprintf(fileID,'Controls: \n');
    fprintf(fileID,'----------------------------------------\n');
    if strcmp(s,'admit')
         fprintf(fileID,'cont.Hmax:  %2i\n',cont.Hmax);
         fprintf(fileID,'cont.IS.opt:  %2i\n',cont.IS.opt);
         fprintf(fileID,'cont.IS.scale:\t');
         fprintf(fileID,'%4.2f\t',cont.IS.scale);
         fprintf(fileID,'\n');
         fprintf(fileID,'cont.IS.perc:\t');
         fprintf(fileID,'%4.2f\t',cont.IS.perc);
         fprintf(fileID,'\n');
         fprintf(fileID,'cont.dfnc:  %2i\n',cont.dfnc);
    else
        fprintf(fileID,'cont.mit.Hmax:  %2i\n',cont.mit.Hmax);
        fprintf(fileID,'cont.mit.dfnc:  %2i\n',cont.mit.dfnc);
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'Preliminary: \n');
    fprintf(fileID,'----------------------------------------\n');
    fprintf(fileID,'Preliminary MH acceptance rate: %6.4f\n',accept);
    fprintf(fileID,'VaR_prelim:  %6.4f\n',VaR_prelim);
    fprintf(fileID,'\n');

    fprintf(fileID,'Mit: \n');
    fprintf(fileID,'----------------------------------------\n');
    if exist('mean_VaR_IS','var')
        fprintf(fileID,'IS (mean) VAR estimate: %6.4f\n',mean_VaR_IS);
    else
        fprintf(fileID,'IS VAR estimate: %6.4f\n',VaR_IS);

    end
    if exist('mean_ES_IS','var')
        fprintf(fileID,'IS (mean) ES estimate: %6.4f\n',mean_ES_IS);
    else
        fprintf(fileID,'IS ES estimate: %6.4f\n',ES_IS);        
    end
    fprintf(fileID,'\n');

    fprintf(fileID,'No. of components q1: %d\n',size(mit1.mu,1));
	fprintf(fileID,'Df. of components q1:\t');
    fprintf(fileID,'%4.2f\t',mit1.df);
    fprintf(fileID,'\n');
	
    fprintf(fileID,'No. of components q2: %d\n',size(mit2.mu,1));
	fprintf(fileID,'Df. of components q2:\t');
    fprintf(fileID,'%4.2f\t',mit2.df);
    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    if exist('summary1','var')
        fprintf(fileID,'CoV q1: %6.4f.\n',summary1.CV(1,end));
    end
    if exist('summary2','var')
        fprintf(fileID,'CoV q2: %6.4f.\n',summary2.CV(1,end));
    end
    fprintf(fileID,'\n');

    fprintf(fileID,'NSE: \n');
    fprintf(fileID,'----------------------------------------\n');
    if exist('NSE_VaR_IS','var')
        fprintf(fileID,'NSE IS VaR estimate: %6.4f.\n',NSE_VaR_IS);
    end
    if exist('NSE_ES_IS','var')
        fprintf(fileID,'NSE IS ES estimate: %6.4f.\n',NSE_ES_IS);
    end
end

fprintf(fileID,'\n');
fclose(fileID);