for ii = 1:11
    figure(ii)
    hist([draw0(:,ii),draw_pmit(:,ii)])
end
kernel = @(x) posterior_debug_hl(x, y, a, b, mean(VaR_prelim), true);
[draw_pmit, lnk_pmit] = fn_p_rmvgt(size(draw0,1), pmit, d, partition, kernel, fn_const_X);


    kernel = @(xx) posterior_debug(xx, y, a, b, true);
    lnk_pmit = kernel(draw_pmit(:,1));

    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    eps_pdf = sum(kernel(draw_pmit(:,2:d)),2);
    lnk_pmit = lnk_pmit + eps_pdf;

lnd_pmit = fn_dpmit(draw_pmit, pmit, partition, fn_const_X, true, GamMat);

w_pmit = lnk_pmit - lnd_pmit;
wmax = max(w_pmit);

w_pmit = w_pmit - wmax;
ew_pmit = exp(w_pmit);



y_hl = predict_arch(draw0(:,1), y_T, S, H, draw0(:,2:end));
PL_hl = sort(fn_PL(y_hl));
  
 boxplot([VaR_prelim,VaR_adapt,VaR_step2,VaR_step2_up],{'VaR prelim','VaR adapt','VaR pmit','VaR pmit up'})
