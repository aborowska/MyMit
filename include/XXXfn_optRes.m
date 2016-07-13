function [mu, Sigma, time] = fn_optIS(theta, w, cont)
    tic

    c = 100*mean(w);
    ISstop = false;
    Sigma = [];
    
    while (IS == false)
        wres = w - c;
        wres(wres<0) = 0;
        [mu, tmp_Sigma] = fn_muSigma(theta, wres); % Sigma - in the vector form
        if (~all(isfinite(tmp_Sigma)) || fn_testSigma(tmp_Sigma))
            c = 0.5*c; % if the matrix not PD or NaNs detected - scale the constant
        else 
            ISstop = true;
        end
    end
    
    for sc = cont.IS.scale
        Sigma = [Sigma; sc*tmp_Sigma];
        % iterate over scaling factors
    end
    time = toc;
end

