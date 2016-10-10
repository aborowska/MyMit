function Print_sv_sml(y,results_path,figures_path,save_on)
    parameter = {'$c$','$\phi$','$\sigma^{2}_{\eta}$','$\nu$'};
 
    if strcmp(results_path, 'results/EQERMit/old/')
        name_add = 'ibm' ;
        time = [2007,2012];
    elseif strcmp(results_path, 'results/EQERMit/crisis/')
        name_add = 'gspc_updated' ;
        time = [2005,2016.5];        
    elseif strcmp(results_path, 'results/EQERMit/recent/')
        name_add = 'gspc_updated_short_sv' ;
        time = [2011 + 2/3,2016 + 2/3];   
    else 
        error('Unknown data.')
    end

    name = ['include/NAIS/SML_',name_add,'.mat'];
    load(name)
    SV_param = par_SV_opt;
    SV_se = sqrt(diag(V_SV_corr_opt));
    SV_theta = theta_smooth;
    
    name = ['include/NAIS/SMLt_',name_add,'.mat'];
    load(name)
    SVt_param = par_SV_opt;
    SVt_se = sqrt(diag(V_SV_corr_opt));
    SVt_theta = theta_smooth;
   
    % Smoothed signal plot
    
    ff = figure(1);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(time(1,1),time(1,2),length(y)); 
    plot(xx,y,'Color',[0.7 0.7 0.7],'LineWidth',1)   
    set(gca,'XLim',time)
    hold on
    plot(xx,SV_theta,'b','LineWidth',2)
    plot(xx,SVt_theta,'Color',[0 139/255 139/255],'LineWidth',2)
    hold off
    
    leg = legend('Data', 'SV', 'SVt');
    set(leg,'Interpreter','latex','FontSize', 12,'Location','NorthEast');
    if save_on
        name = [figures_path,'signal.eps'];
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try                  
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
        end
    end
    
    %% Latex table
    fname = [results_path,'sml.tex']; 
    FID = fopen(fname, 'w+');
    
    fprintf(FID, '\\small{  \n');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.3} \n');
    
    fprintf(FID, '\\begin{longtable}{cc cc c cc } \n');
    caption = ['\\caption{SML estimation results for the parameters of the SV model.} \n'];
    fprintf(FID, caption);

    label = '\\label{tab:sml} \\\\ \n';
    fprintf(FID, label);
    fprintf(FID, '  & &  \\multicolumn{2}{c}{SV} & &  \\multicolumn{2}{c}{SV$t$} \\\\ \\cline{1-1} \\cline{3-4} \\cline{6-7}   \n');
    fprintf(FID, ' Parameter & &  Estimate & St. error &&  Estimate & St. error\\\\ \\cline{1-1} \\cline{3-4} \\cline{6-7}   \n');
   
    for ii = 1:4
        fprintf(FID, '%s & &' ,parameter{ii});
        if (ii < 4)
            fprintf(FID, '%6.4f &  ',SV_param(1,ii));
            fprintf(FID, '$[$ %6.4f $]$ && ',SV_se(ii,1));      
        else
            fprintf(FID, '-- & -- && ');              
        end
        
        fprintf(FID, '%6.4f &  ',SVt_param(1,ii));
        fprintf(FID, '$[$ %6.4f $]$  \\\\ \n',SVt_se(ii,1));
     end

    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\end{longtable} \n');
    fprintf(FID, '} \n');
    fprintf(FID, '} \n');    
    fprintf(FID, '\\normalsize \n');    
    fclose(FID);
end