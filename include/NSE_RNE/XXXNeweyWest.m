function var_nw = NeweyWest(theta,L)
    N = size(theta,1);
    theta = theta - repmat(mean(theta,1),N,1);    
    var_nw = sum(theta.^2,1);
    for ii = 1:L
        w = (L-ii+1)/(L+1);
        gamma =  sum(theta(1+ii:N,:).*theta(1:N-ii,:),1);
        var_nw = var_nw + w*(2*gamma); 
    end
    var_nw = var_nw/(N^2);
end