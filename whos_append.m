matObj = matfile('results/sv_x_VaR_IS_2000_part.mat','Writable',true);
whos('-file','results/sv_x_VaR_IS_2000_part.mat')
matObj.y_opt_h1(:,end+1) = y_opt_h1;
whos('-file','results/sv_x_VaR_IS_2000_part.mat')

 save(['results/sv_x_VaR_IS_2000_part.mat'], 'mit1', 'mit2', 'theta_opt', 'x', 'lnk', 'lnp_T', 'lnd_opt', 'w_opt', 'CV1', 'CV2', 'cont', 'cont2', 'p_bar', 'N', 'M', 'N_sim', 'VaR_prelim','VaR_prelim_MC','ES_prelim', 'VaR_IS', 'ES_IS','-v7.3');
