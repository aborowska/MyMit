function [mu, Sigma, time] = fn_ISopt(theta, w, cont) 
% compute mu and Sigma from the largest weights
% Sigma is a matrix with no. of rows equal to 
% #cont.IS.perc * #cont.IS.scale, default 3*3=9
% each row contains a scale matrix in a vector form, i.e. has d^2 colums
    tic
    [~, ind] = max(w);
    mu0 = theta(ind,:);
    Sigma = [];
    
    for pc = cont.IS.perc % [0.1, 0.15, 0.3] 
        % iterate over percentages of importance weights
        ind = fn_select(w, pc); % get the indeces of draws with the pc% largest IS weights
        [mu, tmp_Sigma] = fn_muSigma(theta(ind,:), w(ind), mu0);
        for sc = cont.IS.scale
            Sigma = [Sigma; sc*tmp_Sigma];
            % iterate over scaling factors
        end
    end
    time = toc;
end