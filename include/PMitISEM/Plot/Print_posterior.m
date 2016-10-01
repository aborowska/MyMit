function results = Print_posterior(results, model, parameter, kernel, GamMat, results_path)
  
    M = 10000;
    BurnIn = 1000;
    N_sim = 20;
    p_bar = 0.01;
    H = 10;

    model_tex = fn_model_tex(model);

    %%
    if isempty(results)
        algo = 'Direct';
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'mit_direct','accept_direct','time_direct')

        theta_mle = mit_direct.mu;
        results.theta_mle = theta_mle;
        d = size(theta_mle,2);
        theta_mle_std = sqrt(diag(reshape(mit_direct.Sigma,d,d)))';
        results.theta_mle_std = theta_mle_std;
        
        theta_direct = Mit_MH(M+BurnIn, kernel, mit_direct, GamMat);
        theta_direct = theta_direct(BurnIn+1:M+BurnIn,:);
        IF_direct = zeros(1,d);
        for ii = 1:d
            [h, ~, bounds] = autocorr(theta_direct(:,ii),100);
            L = min(find((h<bounds(1,1)) & (h>bounds(2,1))));
            if isempty(L)
                L = length(h);
            end
            IF_direct(1,ii) = 1 + 2*sum(h(1:L));
        end
        results.mit_direct = mit_direct;
        results.accept_direct = accept_direct;
%         results.RNE_direct = RNE_direct;
        results.mean_theta_direct = mean(theta_direct);
        results.theta_direct = theta_direct;
        results.IF_direct = IF_direct;
        results.time_direct = time_direct;

        %%
        algo = 'Prelim';
        name = [results_path,model,'_',algo,'_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'mit1','cont1','summary1','accept','time_prelim')
        accept_prelim = accept; 
        clear accept

        theta_prelim = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
        theta_prelim = theta_prelim(BurnIn+1:M+BurnIn,:);

        IF_prelim = zeros(1,d);
        for ii = 1:d
            [h, ~, bounds] = autocorr(theta_prelim(:,ii),100);
            L = min(find((h<bounds(1,1)) & (h>bounds(2,1))));
            if isempty(L)
                L = length(h);
            end
            IF_prelim(1,ii) = 1 + 2*sum(h(1:L));
        end 

        results.mit_prelim = mit1;
        results.accept_prelim = accept_prelim;
%         results.RNE_prelim = RNE_prelim;
        results.mean_theta_prelim = mean(theta_prelim);
        results.theta_prelim = theta_prelim;
        results.IF_prelim = IF_prelim;
        results.time_prelim = time_prelim;
    end
    
    %% Latex table
    fname = [results_path,'results_',model,'_posterior.tex'];

    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');

    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '\\centering \n');

    caption = ['\\caption{Simulation results in the ', model_tex,...
        ' model: estimated posterior means, posterior standard deviations (SD),',...
        ' inefficiency factors (IF) and acceptance rates (AR) of the Markov chains of draws.',...
        ' ML results for comparison.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:posterior_',model,'} \n'];
    fprintf(FID, label);
    fprintf(FID, '\\begin{tabular}{cc rr c rrr c rrr}  \n');
    fprintf(FID, [' & & \\multicolumn{2}{c}{ML}',...
        ' & & \\multicolumn{3}{c}{Naive} & & \\multicolumn{3}{c}{Adapted} ',...
        '\\\\  \\cline{3-4} \\cline{6-8} \\cline{10-12} \n']);   
    fprintf(FID, [' Parameter & &  \\multicolumn{1}{c}{MLE} &  \\multicolumn{1}{c}{SD} ',...
        ' & &  \\multicolumn{1}{c}{Mean} &  \\multicolumn{1}{c}{SD} &  \\multicolumn{1}{c}{IF}',...
        ' & &  \\multicolumn{1}{c}{Mean} &  \\multicolumn{1}{c}{SD} &  \\multicolumn{1}{c}{IF} \\\\ ',...
        '\\cline{1-1}  \\cline{3-4} \\cline{6-8} \\cline{10-12}  \n']);

    for ii = 1:d
        fprintf(FID, '%s & &' ,parameter{ii});

        fprintf(FID,'%7.4f & ', theta_mle(:,ii));
        fprintf(FID,'%7.4f & &', theta_mle_std(:,ii));
        
        fprintf(FID,'%7.4f & ', mean(theta_direct(:,ii)));
        fprintf(FID,'%7.4f & ', std(theta_direct(:,ii)));
        fprintf(FID,'%7.4f & &', IF_direct(:,ii));

        fprintf(FID,'%7.4f & ', mean(theta_prelim(:,ii)));
        fprintf(FID,'%7.4f & ', std(theta_prelim(:,ii)));
        fprintf(FID,'%7.4f \\\\ [1ex] \n', IF_prelim(:,ii));
    end   
    fprintf(FID, '\\cline{1-1}  \\cline{3-4} \\cline{6-8} \\cline{10-12}   \n');

    fprintf(FID, '\\multicolumn{4}{r}{AR} & &');
    fprintf(FID, '\\multicolumn{3}{c}{%6.4f} &&', mean(accept_direct));
    fprintf(FID, '\\multicolumn{3}{c}{%6.4f} \\\\ \n  ', mean(accept_prelim));

    fprintf(FID, '\\cline{1-4} \\cline{6-8} \\cline{10-12}  \n');


    fprintf(FID, ' \\multicolumn{4}{r}{Time construction} & &');
    fprintf(FID, '\\multicolumn{3}{c}{%4.2f s} &&', time_direct(1,1));
    fprintf(FID, '\\multicolumn{3}{c}{%4.2f s} \\\\ \n  ', time_prelim(1,1));

    fprintf(FID, ' \\multicolumn{4}{r}{Time sampling} & &');
    fprintf(FID, '\\multicolumn{3}{c}{%4.2f s} &&', time_direct(2,1));
    fprintf(FID, '\\multicolumn{3}{c}{%4.2f s} \\\\ \n  ', time_prelim(2,1));

    fprintf(FID, ' \\multicolumn{4}{r}{No. of draws }& &');
    fprintf(FID, '\\multicolumn{3}{c}{%s} &&', num2bank(M));
    fprintf(FID, '\\multicolumn{3}{c}{%s} \\\\ \n  ', num2bank(M));
    fprintf(FID, '\\cline{1-4} \\cline{6-8} \\cline{10-12} \n');

    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\end{tabular} \n');
    fprintf(FID, '\\end{table} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '\\normalsize \n');    
    fclose(FID);
end