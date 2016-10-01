function X = WN_ML_const_X(eps, sigma2)
    [N, h] = size(eps);
%     if (h > 1) % not the first forecast
        y = sqrt(sigma2).*eps;
        X = [ones(N,1), fn_PL(y)];
%     else
%         X = ones(N,1);
%     end
end

