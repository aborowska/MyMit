function Plot_hor(y, model, DD, H, p_bar, N_sim, M, BurnIn, save_on, kernel, fn_vol, fn_predict,fn_predict2, partition, fn_const_X, fn_input_X,GamMat)
    close all
    addpath(genpath('include/'));
    d = H + DD;
    y_T = y(end);

    name = ['results/PMitISEM/',model,'_Prelim_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name)
    name = ['results/PMitISEM/',model,'_PMitISEM_',num2str(p_bar),'_H',num2str(H),'_VaR_results_Nsim',num2str(N_sim),'.mat'];
    load(name)

    theta1 = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    theta1 = theta1(BurnIn+1:M+BurnIn,:);   
    vv1 = fn_vol(theta1);
    y_H = fn_predict(theta1,vv1);    
    ind_real = (sum(imag(y_H),2)==0);
    M_real = sum(ind_real); 
    y_H = y_H(ind_real,:);   
    PL_H = sort(fn_PL(y_H));
    VaR_prelim_sim  = PL_H(round(p_bar*M_real));
    Plot_hor_direct(y_H, y_T, VaR_prelim_sim, model, save_on)


    theta_pmit  = fn_p_rmvgt2(M/2, pmit, d, partition, [], fn_const_X, fn_input_X);      
%     f_T = volatility_t_gas_mex(draw_pmit(:,1:DD), y);
%     y_pmit  = predict_t_gas(draw_pmit(:,1:DD), y_T, f_T, H, draw_pmit(:,DD+1:H+DD));
    vv_pmit = fn_vol(theta_pmit(:,1:DD));
    y_pmit = fn_predict2(theta_pmit(:,1:DD),vv_pmit,theta_pmit(:,DD+1:H+DD));
    Plot_hor_pmit(y_pmit,y_T, VaR_prelim, model,'PMitISEM', save_on)
end