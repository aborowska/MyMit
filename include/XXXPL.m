function [val_return] = fn_PL(vars, L)
% profit loss function value at c 
% if natin>1
% if L == 1 ==> VaR 
% if L == 2 ==> PL density estimator at c
% w - weights corresponding to y

    f_pl = @(y) 100*(exp(y/100) - 1);
    
    if (nargin > 1) % if there is L
        y = vars.y;
        w = vars.w;
       
        PL = f_pl(y);
        [PL, ind] = sort(PL); 
        w = w(ind,:);
        w = w/sum(w);
        cum_w = cumsum(w);
        
        if (L == 1) % compute VaR_IS and ES_IS
            % i.e.: find c s.t. p_hat(PL(X)<=c)=p_bar 
            p_bar = vars.p_bar;
            ind_var = min(find(cum_w > p_bar))-1; 
            VaR_IS = PL(ind_var);
            ES_IS = sum((w(1:ind_var)/sum(w(1:ind_var))).*PL(1:ind_var));
            val_return = [VaR_IS, ES_IS]; 
        elseif (L == 2) % compute p_hat(c)
            c = vars.c;
            p_hat = zeros(10,1);
            % c == VaR_IS
            for eps=0.01:0.01:0.1
                ind_p = min(find(PL > c + eps))-1;
                ind_m = min(find(PL > c - eps))-1;
                p_hat(floor(eps*100),1) = (cum_w(ind_p) - cum_w(ind_m))./(2*eps);
            end
            p_hat = min(p_hat(p_hat~=0));
            val_return = p_hat;
        end
    elseif (nargin == 1)
        c = vars;
        pl = f_pl(c);
        val_return = pl;
    end
    
    
end

