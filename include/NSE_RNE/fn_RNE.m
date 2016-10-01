function RNE = fn_RNE(theta, algo, input, weighting,p)
% Computes Relative Numerical Efficiency for the simulated sample theta
% theta can come from the MH algorithm (algo: 'MH') or from importance sampling (algo:'IS')
% For 'MH' input is an optional lag parameter (default: the rule of thumb with L= 0.75*N^(1/3))
% For 'IS' input are importance weights (optimal, default: equal)
% For 'MH' three weighting schemes:  standard Newey-West weights ('NW'), 
% Parzen kernel ('P') or (default) Quadratic kernel ('Q') 
% For 'MH' performs prewhitening with the default filter order p=1
    [N,d] = size(theta);
    
    if strcmp(algo,'MH')
        %  the quality of a correlated sample
        % It compares the empirical variance of the sample with a correlation-consistent t variance estimator
        if ((nargin == 2) || isempty(input))
            input = floor(0.75*N^(1/3)-1);
        end
        if ((nargin < 4) || isempty(weighting))
            weighting = 'Q';
        end
        if (nargin < 5) 
            p = 1;
        end
        % Newey-West vs. iid
        var_direct = var(theta)/N; 
        var_nw = NeweyWest(theta,input,weighting,p); 
        % This estimator tends to underestimate the autocovariance matrices
        % especially for larger values of aurocorr order
        
        RNE = var_direct/var_nw;
    elseif strcmp(algo,'IS')
        % Direct sampling 
        if (nargin == 3)
            w = input;
        else
            w = ones(N,1);
        end
        w = w/sum(w);
        tmp_w = repmat(w,1,d);
        theta_mean = sum(tmp_w.*theta,1);
        theta_demeaned = theta - repmat(theta_mean,N,1);
        var_direct = (sum(tmp_w.*(theta.^2),1) - theta_mean.^2)/N;
        var_is = sum((tmp_w.^2).*(theta_demeaned.^2),1); %NSE = sqrt(var_is);
        RNE = var_direct/var_is;
    end
end

function var_nw = NeweyWest(theta,L,weighting,p)
    theta = PreWhiten(theta,p);
    N = size(theta,1);
    
    if strcmp(weighting,'P') % parzen kernel
        PR_weights = @(aaa) (1 - 6*aaa.^2 + 6*abs(aaa).^3).*(abs(aaa)<=0.5) ...
                    +  (2*(1-abs(aaa)).^3).*(abs(aaa)>0.5).*(abs(aaa)<=1);
        w = PR_weights((1:L)/L);
    elseif strcmp(weighting,'Q') % quadratic spectral
        QS_weights = @(aaa) 25*(sin(6*pi*aaa/5)./(6*pi*aaa/5) - cos(6*pi*aaa/5))./(12*(pi^2).*(aaa.^2));
        w = QS_weights((1:L)/L);
    else % standard NW weights
        w = (L+1-(1:L))/(L+1);
    end
    
    theta = theta - repmat(mean(theta,1),N,1);    
    var_nw = sum(theta.^2,1);
    for ii = 1:L
        gamma =  sum(theta(1+ii:N,:).*theta(1:N-ii,:),1);
        var_nw = var_nw + w(1,ii)*(2*gamma); 
    end
    var_nw = var_nw/(N^2);
end


function theta_pw = PreWhiten(theta,p)
    if p > 0
        [N,d] = size(theta);
        theta_pw = zeros(N-p,d);    
        Y = theta(p+1:N,:);
        for ii = 1:d
            X = zeros(N-p,p);
            for jj = 1:p
                X(:,jj) = theta(p+1-jj:N-jj,ii);
            end
            beta = (X'*X)\X'*Y;
            theta_pw(:,ii) = Y - X*beta;
        end
    else
        theta_pw = theta;
    end
end