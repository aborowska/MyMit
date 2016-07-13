fileID = fopen(['results/',model,'_',algo,'_',num2str(p_bar),'_',num2str(M),'.txt'],'w');

fprintf(fileID,'Model: %s \n',model);
fprintf(fileID,'Algorithm: %s \n',algo);
fprintf(fileID,'----------------------------------------\n');

fprintf(fileID,'Controls: \n');
fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'N_sim:  %2i\n',N_sim);
if strcmp(algo,'MitISEM')
    fprintf(fileID,'N:  %2i\n',cont.mit.N);
elseif strcmp(algo,'AdMit')
    fprintf(fileID,'Ns:  %2i\n',cont.Ns);
    fprintf(fileID,'Np:  %2i\n',cont.Np);
end
fprintf(fileID,'M:  %2i\n',M);
if strcmp(algo,'MitISEM')
    fprintf(fileID,'df1:  %2i\n',cont.mit.dfnc);
    fprintf(fileID,'df2:  %2i\n',cont2.mit.dfnc);
elseif strcmp(algo,'AdMit')
    fprintf(fileID,'df1:  %2i\n',cont.dfnc);
    fprintf(fileID,'df2:  %2i\n',cont2.dfnc);
else
    fprintf(fileID,'df:  %2i\n',cont.mit.dfnc);
end
fprintf(fileID,'\n');
fprintf(fileID,'Acceptance: \n');
fprintf(fileID,'----------------------------------------\n');
if strcmp(algo,'Direct')
    fprintf(fileID,'mean MH acceptance:  %6.4f.\n',mean(accept_direct));
else 
    fprintf(fileID,'mean MH acceptance:  %6.4f.\n',mean(accept));
end
fprintf(fileID,'\n');
fprintf(fileID,'Estimates: \n');
fprintf(fileID,'----------------------------------------\n');
if ~strcmp(algo,'Direct')
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
else
    fprintf(fileID,'100*%4.2f%% VaR direct (mean) estimate: %6.4f. \n', p_bar, mean_VaR_direct);
    fprintf(fileID,'NSE VaR direct: %6.4f. \n',  NSE_VaR_direct);
    fprintf(fileID,'VaR direct: [%6.4f, %6.4f]. \n \n',  mean_VaR_direct - NSE_VaR_direct, mean_VaR_direct + NSE_VaR_direct);


    fprintf(fileID,'100*%4.2f%% ES direct (mean) estimate: %6.4f. \n',  p_bar, mean_ES_direct);
    fprintf(fileID,'NSE ES direct: %6.4f. \n',  NSE_ES_direct);
    fprintf(fileID,'ES direct: [%6.4f, %6.4f]. \n \n', mean_ES_direct - NSE_ES_direct, mean_ES_direct + NSE_ES_direct);
end
fclose(fileID);
