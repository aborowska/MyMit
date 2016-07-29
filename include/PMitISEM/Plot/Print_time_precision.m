function results = Print_time_precision(model,horizons,p_bar,N_sim,M,estimation)
    
    if isempty(strfind(model,'_ML'))
        ML = false;
    else
        ML = true;
    end
    
    switch model
        case 't_gas'
            model_tex = 'GAS(1,1)-$t$';  
        case 't_gas_ML'
            model_tex = 'GAS(1,1)-$t$';        
        case 't_garch2_noS'
            model_tex = 'GARCH(1,1)-$t$';
        case 'arch'
            model_tex = 'ARCH(1)';
        case 'WN'
            model_tex = 'White Noise';
        case 'WN_ML'
            model_tex = 'White Noise';
    end
  
    if (nargin < 6)
        estimation = '';
    else
        model_tex = ['WN(',estimation,')'];
        estimation = [estimation,'_'];
    end
    
    Hno = size(horizons,2);
 
    VaR_NSE = zeros(Hno,4);
    ES_NSE = zeros(Hno,4);
    VaR_precision = zeros(Hno,4);
    ES_precision = zeros(Hno,4);
    time_construction = zeros(Hno,4);
    time_sampling = zeros(Hno,4);
    h = 0;
    
    for H = horizons
        h = h+1;
        % load results
        name = ['results/PMitISEM/',model,'_Direct_',estimation,num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        if ~ML
          	load(name,'VaR_direct','ES_direct','time_direct')
        else
            load(name,'VaR_direct','ES_direct','time_direct','time_bigdraw')
        end
        VaR_NSE(h,1) = std(VaR_direct);
        ES_NSE(h,1) = std(ES_direct);
        VaR_precision(h,1) = 1/var(VaR_direct);
        ES_precision(h,1) = 1/var(ES_direct);        
        time_construction(h,1) = time_direct(1,1);
        time_sampling(h,1) = time_direct(2,1);
        clear VaR_direct ES_direct time_direct
        
        if ~ML
            name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
            load(name,'VaR_prelim','ES_prelim','time_prelim','time_bigdraw') 
            prelim_ok = true;
        else
            VaR_prelim = NaN*ones(N_sim,1);
            ES_prelim = NaN*ones(N_sim,1);
            time_prelim = zeros(2,1);
            prelim_ok = false;
        end
        VaR_NSE(h,2) = std(VaR_prelim);
        ES_NSE(h,2) = std(ES_prelim);
        VaR_precision(h,2) = 1/var(VaR_prelim);
        ES_precision(h,2) = 1/var(ES_prelim);        
        time_construction(h,2) = time_prelim(1,1);
        time_sampling(h,2) = time_prelim(2,1);
        clear VaR_prelim ES_prelim time_prelim        
        
        try
            name = ['results/PMitISEM/',model,'_MitISEM_',estimation,num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
            load(name,'VaR_mit','ES_mit','time_mit') 
            mit_ok = true;    
            if ~exist('VaR_mit','var')
                VaR_mit = NaN*ones(N_sim,1);
                ES_mit = NaN*ones(N_sim,1);
                time_mit = NaN*ones(2,1);
                mit_ok = false;
            end
        catch
            VaR_mit = NaN*ones(N_sim,1);
            ES_mit = NaN*ones(N_sim,1);
            time_mit = NaN*ones(2,1);
            mit_ok = false;
        end

        VaR_NSE(h,3) = std(VaR_mit);
        ES_NSE(h,3) = std(ES_mit);
        VaR_precision(h,3) = 1/var(VaR_mit);
        ES_precision(h,3) = 1/var(ES_mit);        
        time_construction(h,3) = time_mit(1,1) + time_construction(h,2);
        time_sampling(h,3) = time_mit(2,1);
        clear VaR_mit ES_mit time_mit
        
        name = ['results/PMitISEM/',model,'_PMitISEM_',estimation,num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
        load(name,'VaR_pmit','ES_pmit','time_pmit')    
        VaR_NSE(h,4) = std(VaR_pmit);
        ES_NSE(h,4) = std(ES_pmit);
        VaR_precision(h,4) = 1/var(VaR_pmit);
        ES_precision(h,4) = 1/var(ES_pmit);        
        time_construction(h,4) = time_pmit(1,1) + time_construction(h,2);
        time_sampling(h,4) = time_pmit(2,1);
        clear VaR_pmit ES_pmit time_pmit     
    end
    
    if ML
        VaR_NSE(:,2) = [];
        ES_NSE(:,2) = [];
        VaR_precision(:,2) = [];
        ES_precision(:,2) = [];
        time_construction(:,2) = [];
        time_sampling(:,2) = [];
    end
    precision_one_digit = (1.96/0.05)^2; % 1.96*NSE<=0.05

    VaR_slope = VaR_precision./time_sampling;
    ES_slope = ES_precision./time_sampling;

    VaR_intercept = -VaR_slope.*time_construction;
    ES_intercept = -ES_slope.*time_construction;

    time_total = time_construction + time_sampling;
    time_per_draw = M./time_sampling;
        
    VaR_time_required = (precision_one_digit - VaR_intercept)./VaR_slope;   
    ES_time_required = (precision_one_digit - ES_intercept)./ES_slope;   
       
    VaR_draws_required = M*precision_one_digit./VaR_precision;   
    ES_draws_required = M*precision_one_digit./ES_precision;   
    
    
    
    results.time_construction = time_construction;
    results.time_sampling = time_sampling;
    results.time_total = time_total;
    results.time_per_draw = time_per_draw; 
    
    results.VaR_NSE = VaR_NSE;
    results.ES_NSE = ES_NSE;
    
    results.VaR_precision = VaR_precision;
    results.ES_precision = ES_precision;
    
    results.VaR_slope = VaR_slope;
    results.ES_slope = ES_slope;
    
    results.VaR_intercept = VaR_intercept; 
    results.ES_intercept = ES_intercept;
        
    results.VaR_time_required = VaR_time_required; 
    results.ES_time_required = ES_time_required;
       
    results.VaR_draws_required = VaR_draws_required;  
    results.ES_draws_required = ES_draws_required; 
 
    
    %% create a tex table
    fname = ['results/PMitISEM/results_',model,'_',estimation,'time_precision_comb.tex'];
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '\\centering \n');

    caption = ['\\caption{Computation time-precision trade-off for the  $99\\%%$ VaR and ES evaluation in ', model_tex,' model for different horizons.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:time_precision_',model,'} \n'];
    fprintf(FID, label);
    if ML 
        fprintf(FID, '\\begin{tabular}{rr rrr r rrr}  \n');
        fprintf(FID, ' & & \\multicolumn{1}{c}{Direct} & \\multicolumn{2}{c}{QERMit}&  & \\multicolumn{1}{c}{Direct} & \\multicolumn{2}{c}{QERMit} \\\\ \\cline{3-5} \\cline{7-9} \n');   
        fprintf(FID, ' H & & Naive & MitISEM & PMitISEM & & Naive & MitISEM & PMitISEM \\\\ \\hline \n');       
    else
        fprintf(FID, '\\begin{tabular}{rr rrrr r rrrr}  \n');
        fprintf(FID, ' & & \\multicolumn{2}{c}{Direct} & \\multicolumn{2}{c}{QERMit}&  & \\multicolumn{2}{c}{Direct} & \\multicolumn{2}{c}{QERMit} \\\\ \\cline{3-6} \\cline{8-11} \n');   
        fprintf(FID, ' H & & Naive & Adapted & MitISEM & PMitISEM & & Naive & Adapted & MitISEM & PMitISEM \\\\ \\hline \n');
    end
    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{Total time}  \\\\ \\cline{3-5} \n');     
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{Total time}  \\\\ \\cline{3-6} \n');     
    end
    for h = 1:Hno
        if ML        
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s  \\\\ \n',horizons(1,h), time_total(h,:));            
        else
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n',horizons(1,h), time_total(h,:));
        end
    end
    fprintf(FID, '\\hline \n');    
    
    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{Construction time} & & \\multicolumn{3}{c}{ Sampling time} \\\\ \\cline{3-5}  \\cline{7-9}\n');         
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{Construction time} & & \\multicolumn{4}{c}{ Sampling time} \\\\ \\cline{3-6}  \\cline{8-11}\n');
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s &&  %4.2f s & %4.2f s & %4.2f s \\\\ \n',horizons(1,h), time_construction(h,:), time_sampling(h,:));            
        else
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s & %4.2f s && %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n',horizons(1,h), time_construction(h,:), time_sampling(h,:));
        end
    end
    fprintf(FID, '\\hline \n');

    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{VaR NSE} &&  \\multicolumn{3}{c}{ES NSE} \\\\ \\cline{3-5}  \\cline{7-9}\n');             
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{VaR NSE} &&  \\multicolumn{4}{c}{ES NSE} \\\\ \\cline{3-6}  \\cline{8-11}\n'); 
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i && %6.4f  & %6.4f  & %6.4f && %6.4f  & %6.4f  & %6.4f  \\\\ \n',horizons(1,h), VaR_NSE(h,:),  ES_NSE(h,:));
        else
            fprintf(FID, '%i && %6.4f  & %6.4f  & %6.4f & %6.4f && %6.4f  & %6.4f  & %6.4f & %6.4f \\\\ \n',horizons(1,h), VaR_NSE(h,:),  ES_NSE(h,:));
        end
    end
    fprintf(FID, '\\hline \n');

    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{VaR precision} &&  \\multicolumn{3}{c}{ES precision} \\\\ \\cline{3-5}  \\cline{7-9}\n');             
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{VaR precision} &&  \\multicolumn{4}{c}{ES precision} \\\\ \\cline{3-6}  \\cline{8-11}\n');     
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i &&  %6.4f & %6.4f & %6.4f & & %6.4f & %6.4f & %6.4f \\\\ \n',horizons(1,h), VaR_precision(h,:),  ES_precision(h,:));
        else            
            fprintf(FID, '%i && %6.4f & %6.4f & %6.4f & %6.4f & & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n',horizons(1,h), VaR_precision(h,:),  ES_precision(h,:));
        end
    end
    fprintf(FID, '\\hline \n');

    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{ VaR slope} && \\multicolumn{3}{c}{ES slope} \\\\ \\cline{3-5}  \\cline{7-9}\n');         
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{ VaR slope} && \\multicolumn{4}{c}{ES slope} \\\\ \\cline{3-6}  \\cline{8-11}\n'); 
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i && %4.2f & %4.2f & %4.2f && %4.2f & %4.2f & %4.2f \\\\ \n',horizons(1,h), VaR_slope(h,:),  ES_slope(h,:));            
        else
            fprintf(FID, '%i && %4.2f & %4.2f & %4.2f & %4.2f && %4.2f & %4.2f & %4.2f & %4.2f \\\\ \n',horizons(1,h), VaR_slope(h,:),  ES_slope(h,:));
        end
    end
    fprintf(FID, '\\hline \n');

    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{ VaR intercept} &&  \\multicolumn{3}{c}{ES intercept} \\\\ \\cline{3-5}  \\cline{7-9}\n');     
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{ VaR intercept} &&  \\multicolumn{4}{c}{ES intercept} \\\\ \\cline{3-6}  \\cline{8-11}\n');
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i &&  %4.2f & %4.2f & %4.2f && %4.2f & %4.2f & %4.2f \\\\ \n',horizons(1,h), VaR_intercept(h,:),  ES_intercept(h,:));
        else            
            fprintf(FID, '%i && %4.2f & %4.2f & %4.2f & %4.2f && %4.2f & %4.2f & %4.2f & %4.2f \\\\ \n',horizons(1,h), VaR_intercept(h,:),  ES_intercept(h,:));
        end
    end
    fprintf(FID, '\\hline \n');

    
    if ML
        fprintf(FID, ' & & \\multicolumn{3}{c}{VaR time required} && \\multicolumn{3}{c}{ES time required} \\\\ \\cline{3-5}  \\cline{7-9}\n');     
    else
        fprintf(FID, ' & & \\multicolumn{4}{c}{VaR time required} && \\multicolumn{4}{c}{ES time required} \\\\ \\cline{3-6}  \\cline{8-11}\n');     
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s && %4.2f s & %4.2f s & %4.2f s \\\\ \n',horizons(1,h), VaR_time_required(h,:), ES_time_required(h,:));
        else
            fprintf(FID, '%i & & %4.2f s & %4.2f s & %4.2f s & %4.2f s && %4.2f s & %4.2f s & %4.2f s & %4.2f s \\\\ \n',horizons(1,h), VaR_time_required(h,:), ES_time_required(h,:));
        end
    end
    fprintf(FID, '\\hline \n');
    
    if ML
        fprintf(FID, ' && \\multicolumn{3}{c}{VaR draws required} &&   \\multicolumn{3}{c}{ES draws required} \\\\  \\cline{3-5}  \\cline{7-9} \n');     
    else
        fprintf(FID, ' && \\multicolumn{4}{c}{VaR draws required} &&   \\multicolumn{4}{c}{ES draws required} \\\\  \\cline{3-6}  \\cline{8-11} \n');     
    end
    for h = 1:Hno
        if ML
            fprintf(FID, '%i & & %i & %i & & %i & %i & %i & %i \\\\ \n',horizons(1,h),  round(VaR_draws_required(h,:)), round(ES_draws_required(h,:)));
        else
            fprintf(FID, '%i & & %i & %i & %i & & %i & %i & %i & %i  & %i \\\\ \n',horizons(1,h),  round(VaR_draws_required(h,:)), round(ES_draws_required(h,:)));
        end
    end
    fprintf(FID, '\\hline \n');
  
        
    fprintf(FID, '\\end{tabular} \n');
    fprintf(FID, '\\end{table} \n');
    fprintf(FID, '} \n');
    fclose(FID);
    
end