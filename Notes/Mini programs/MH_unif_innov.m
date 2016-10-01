function [x_chain, AR] = MH_unif_innov(N,ll,innov)
    u = rand(N,1);
    x_chain = zeros(N,1);
    x_chain(1,1) = innov(1,1);
    old_ll = ll(x_chain(1,1));

    AR = 0;
    for ii = 2:N
        x_cand = x_chain(ii-1,1) + innov(ii,1);
        new_ll = ll(x_cand);
        A = min(1,exp(new_ll - old_ll));
        if (A > u(ii,1))
            x_chain(ii,1) = x_cand;
            AR = AR + 1; 
        else
            x_chain(ii,1) = x_chain(ii-1,1);    
        end
    end
    AR = AR/(N-1);
end