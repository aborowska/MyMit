if strcmp(model,'WN')
    fileID = fopen(['results/',model,'_',num2str(p_bar),'_a_',num2str(a),'_',num2str(M),'.txt'],'w');
else
    fileID = fopen(['results/',model,'_',num2str(p_bar),'_',num2str(M),'.txt'],'w');
end

fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'Model %s: \n',model);
fprintf(fileID,'----------------------------------------\n');

fprintf(fileID,'Controls: \n');
fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'N_sim1:  %2i\n',N_sim1);
fprintf(fileID,'N_sim:  %2i\n',N_sim);
fprintf(fileID,'N:  %2i\n',cont.mit.N);
fprintf(fileID,'M:  %2i\n',M);
fprintf(fileID,'df1:  %2i\n',cont.mit.dfnc);
fprintf(fileID,'df2:  %2i\n',cont2.mit.dfnc);
fprintf(fileID,'\n');
fprintf(fileID,'Estimates: \n');
fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'100*%4.2f%% VaR prelim (mean) estimate: %6.4f. \n', p_bar, mean_VaR_prelim);
fprintf(fileID,'NSE VaR prelim: %6.4f. \n', NSE_VaR_prelim);
fprintf(fileID,'VaR prelim: [%6.4f, %6.4f]. \n \n', mean_VaR_prelim - NSE_VaR_prelim, mean_VaR_prelim + NSE_VaR_prelim);

fprintf(fileID,'100*%4.2f%% VaR IS (mean) estimate: %6.4f. \n', p_bar, mean_VaR_IS);
fprintf(fileID,'NSE VaR IS estimate: %6.4f. \n', NSE_VaR_IS);
fprintf(fileID,'VaR: [%6.4f, %6.4f]. \n \n', mean_VaR_IS - NSE_VaR_IS, mean_VaR_IS + NSE_VaR_IS);

fprintf(fileID,'100*%4.2f%% ES prelim (mean) estimate: %6.4f. \n', p_bar, mean_ES_prelim);
fprintf(fileID,'NSE ES prelim: %6.4f. \n', NSE_ES_prelim);
fprintf(fileID,'ES prelim: [%6.4f, %6.4f]. \n \n', mean_ES_prelim - NSE_ES_prelim, mean_ES_prelim + NSE_ES_prelim);

fprintf(fileID,'100*%4.2f%% ES IS (mean) estimate: %6.4f. \n', p_bar, mean_ES_IS);
fprintf(fileID,'NSE ES IS estimate: %6.4f. \n', NSE_ES_IS);
fprintf(fileID,'ES: [%6.4f, %6.4f]. \n',  mean_ES_IS - NSE_ES_IS, mean_ES_IS + NSE_ES_IS);

fclose(fileID);
