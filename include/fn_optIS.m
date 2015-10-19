function [mu, Sigma] = fn_optIS(theta, w, cont, mu0, mu_mat) 
% compute mu and Sigma from the largest weights
% Sigma is a matrix with no. of rows equal to 
% #cont.IS.perc * #cont.IS.scale, default 3*3=9
% each row contains a scale matrix in a vector form, i.e. has d^2 colums

% mu_mat === true if mu has to be in a matrix form, with rows corresponding
% to different percentages of the higherst weights

    if (nargin == 3)
        [~, ind] = max(w);
        mu0 = theta(ind,:);
    end
    
    if ((nargin == 5) && mu_mat)
        mu = zeros(size(cont.IS.perc,2)*size(cont.IS.scale,2),size(theta,2));    
    end
    Sigma = zeros(size(cont.IS.perc,2)*size(cont.IS.scale,2),size(theta,2)^2);
    ii = 0;
    
    for pc = cont.IS.perc % [0.1, 0.15, 0.3] 
        % iterate over percentages of importance weights
        ind = fn_select(w, pc); % get the indeces of draws with the pc% largest IS weights
        if ((nargin == 5) && mu_mat)
            [tmp_mu, tmp_Sigma] = fn_muSigma(theta(ind,:), w(ind));            
        else
            [mu, tmp_Sigma] = fn_muSigma(theta(ind,:), w(ind), mu0);
        end
        for sc = cont.IS.scale % [1, 0.25, 4]
%             Sigma = [Sigma; sc*tmp_Sigma];
            ii = ii + 1;
            if ((nargin == 5) && mu_mat)
                mu(ii,:) = tmp_mu;
            end
            Sigma(ii,:) = sc*tmp_Sigma;
            % iterate over scaling factors
        end
    end
end