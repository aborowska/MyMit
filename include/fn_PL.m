function [val_return] = fn_PL(vars, L)
% profit loss function value at c 
% if nargin>1
% if L == 1 ==> VaR esimtation via IS
% if L == 2 ==> PL density estimator at c (to explicitly compute NSE of VaR)

% w - weights corresponding to y

    f_pl = @(y) 100*(exp(y/100) - 1); % percentage profit loss
    
    if (nargin > 1) % if there is L
        y = vars.y;
        w = vars.w;
       
        PL = f_pl(sum(y,2));
        [PL, ind] = sort(PL); 
        w = w(ind,:);
        w = w/sum(w);
        cum_w = cumsum(w);
        
        if (L == 1) % compute VaR_IS and ES_IS
            % i.e.: find c s.t. p_hat(PL(X)<=c)=p_bar 
            p_bar = vars.p_bar;
            ind_var = min(find(cum_w > p_bar))-1; 
            if isempty(ind_var)
                ind_var = length(y)-1;
            end
            if (ind_var == 0)
                ind_var = 1;
            end
            VaR_IS = (PL(ind_var+1) + PL(ind_var))/2; % intrapolate
            ES_IS = sum((w(1:ind_var)/sum(w(1:ind_var))).*PL(1:ind_var));
            val_return = [VaR_IS, ES_IS]; 
        elseif (L == 2) % compute p_hat(c)
            c = vars.c;
%             scale = vars.scale;
            scale = 1;
            val_return = [];
            iter_max = 5;
            iter = 0;
            while ((isempty(val_return)) && (iter <= iter_max))
                iter = iter + 1;
                scale = scale/(10^(iter-1));
                p_hat = zeros(10,1);
                % c == VaR_IS                
                for ii = 1:10
                    eps = scale*ii;
                    ind_p = min(find(PL > c + eps))-1;
                    ind_m = max(1,min(find(PL > c - eps))-1);
                    p_hat(ii,1) = (cum_w(ind_p) - cum_w(ind_m))./(2*eps);
                end
                p_hat = min(p_hat(p_hat~=0));
                val_return = p_hat;
            end
        end
    elseif (nargin == 1)
        c = vars;
        val_return = f_pl(sum(c,2));
    end
    
    
end

