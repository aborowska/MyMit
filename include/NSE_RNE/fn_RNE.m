function RNE = fn_RNE(input, algo)
    theta = input.theta;
    [N,d] = size(theta);
    if strcmp(algo,'MH')
        %  the quality of a correlated sample
        % It compares the empirical variance of the sample with a correlation-consisten t variance estimator
        lag = input.lag;
        % Newey-West vs. iid
        var_direct = var(theta);
        var_nw = NeweyWest(theta,N,lag);
        RNE = var_direct/var_nw;
    elseif strcmp(algo,'IS')
        % Direct sampling 
        w = input.w;
        w = w/sum(w);
        tmp_w = repat(w,1,d);
        theta_mean = sum(tmp_w.*theta,1);
        theta_demeaned = theta - repmat(theta_mean,N,1);
        var_direct = (sum(tmp_w.*(theta.^2),1) - theta_mean.^2)/N;
        var_is = sum((tmp_w.^2).*(theta_demeaned.^2),1); %NSE = sqrt(var_is);
        RNE = var_direct/var_is;
    end
end

function var_nw = NeweyWest(theta,N,L)
    theta = theta - repmat(mean(theta,1),N,1);    
    var_nw = theta'*theta;
    for ii = 1:L
        w = (L-ii+1)/(L+1);
        gamma =  theta(1+ii:N,:)'*theta(1:N-ii,:);
        var_nw = var_nw + w*(gamma + gamma'); 
    end
end
 