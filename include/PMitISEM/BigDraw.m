function [draw_mm, VaR, y_H, PL] = BigDraw(M, H, BurnIn, p_bar, mit1, kernel, y_predict, GamMat, nu_col)

    F = 1/p_bar;
    % M = 10000;
    MM = 2*M;
    % kernel = @(x) posterior_debug(x, y, a, b, true);
    fprintf('\nBigDraw\n')
    theta1 = Mit_MH(F*M+BurnIn, kernel, mit1, GamMat);
    theta1 = theta1((BurnIn+1):(F*M+BurnIn),:);
    
	theta_new = theta1(1:MM,:);
    if (nargin == 9) % the optional last argument says in which culumn theta are degrees of freedom
        nu = theta_new(:,nu_col);
        eps_H = trnd(repmat(nu,1,H));
    else
        eps_H = randn(MM,H); 
    end

    draw_mm = [theta_new,eps_H];
    % y_H = bsxfun(@times,eps_H,sqrt(theta_new)); 
    y_H = y_predict(draw_mm); 

    [~, ind] = sort(fn_PL(y_H)); 
    ind_good = ind(1:M);
    ind_bad = ind(M+1:MM);
    % ind_good_old = ind_good(ind_good<=M);
    ind_good_new = ind_good(ind_good>M);
    ind_bad_old = ind_bad(ind_bad<=M);
    draw_mm(ind_bad_old,:) = draw_mm(ind_good_new,:);


    for ii = 3:F
        fprintf('\nBigDraw iter: %s\n',num2str(ii));
        ind_mm = (1:M)' + (ii-1)*M; 
        theta_new = theta1(ind_mm,:);
        
        if (nargin == 9) % the optional last argument says in which culumn theta are degrees of freedom
            nu = theta_new(:,nu_col);
            eps_H = trnd(repmat(nu,1,H));
        else
            eps_H = randn(M,H); 
        end
    %     eps_H = randn(M,H); 
        draw_mm(M+1:MM,:) = [theta_new,eps_H];

%         y_H = bsxfun(@times,draw_mm(:,2:H+1),sqrt(draw_mm(:,1))); 
        y_H = y_predict(draw_mm); 
        [~, ind] = sort(fn_PL(y_H)); 
        ind_good = ind(1:M);
        ind_bad = ind(M+1:MM);
        % ind_good_old = ind_good(ind_good<=M);
        ind_good_new = ind_good(ind_good>M);
        ind_bad_old = ind_bad(ind_bad<=M);
%         fprintf('\n Length ind_good_new: %d\n',length(ind_good_new))
        draw_mm(ind_bad_old,:) = draw_mm(ind_good_new,:);
    end

    % y_H = bsxfun(@times,draw_mm(:,2:H+1),sqrt(draw_mm(:,1)));
    y_H = y_predict(draw_mm); 
    [PL, ind] = sort(fn_PL(y_H)); 
    VaR = PL(M);
    draw_mm = draw_mm(ind,:);
    draw_mm = draw_mm(1:M,:);
    y_H = y_H(ind);
    y_H = y_H(1:M);
    PL = PL(1:M);
end