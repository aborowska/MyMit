function var_nw = NeweyWest(theta,L)
    N = size(theta,1);
    theta = theta - repmat(mean(theta,1),N,1);    
    var_nw = theta'*theta;
    for ii = 1:L
        w = (L-ii+1)/(L+1);
        gamma =  theta(1+ii:N,:)'*theta(1:N-ii,:);
        var_nw = var_nw + w*(gamma + gamma'); 
    end
end