function r = fn_testSigma(vS)
    d = sqrt(length(vS));
    r = fn_isPDS(reshape(vS,d,d));
end

function r = fn_isPDS(M)
% r = true if any problem with PSD of M
% r = false if no problems
    r = any(eig(M)<=0); % if true, M not pd matrix
    if (~r) % if M pd
        try
            r = chol(M);
            r = false;
        catch
            r = true; % problem with Cholesky decomposition
        end
    end
    tmp = abs(det(M));
    tol = 1e25;
    
    if (~r)
        r = ((tmp >= tol) || (tmp <= 1/tol)); % determinant too small/large
    end 
    
    if (~r)
        try
            r = inv(M);
            r = false;
        catch
            r = true; % problem with inversion
        end        
    end
end