function Mit = MitISEM_plot(mit, D, x, y, GamMat)
    [H,d] = size(mit.mu);
 
    if (D == 2)  % 2D plot   
        x = x';
        n = length(x);
        Mit = zeros(n,1);
        for h = 1:H
            p_h = mit.p(h);
            mu_h = mit.mu(h,:);
            Sigma_h = mit.Sigma(h,:);
            df_h = mit.df(h);

            MitX_h = @(a) p_h*dmvt(a, mu_h, Sigma_h, df_h, GamMat);
            MitX_h = arrayfun(@(ii) MitX_h(x(ii,:)), 1:n, 'un', 0);
            MitX_h = cell2mat(MitX_h);
            MitX_h = MitX_h';
            Mit = Mit + MitX_h;
        end
        plot(x,Mit)
        
    elseif (D == 3) % 3D plot
        
        n = length(x);
        m = length(y);
        [X1,X2] = meshgrid(x,y);
        V1 = reshape(X1,n*m,1); V2 = reshape(X2,n*m,1);
        V = [V1,V2];
        Mit = zeros(m,n);

        for h=1:H
            p_h = mit.p(h);
            mu_h = mit.mu(h,:);
            Sigma_h = reshape(mit.Sigma(h,:),d,d);
            df_h = mit.df(h);

            MitX_h = @(a) p_h*dmvt(a, mu_h, Sigma_h, df_h, GamMat);
        %     MitX_h = @(a) 1*dmvt(a,mit_init.mu,reshape(mit_init.Sigma,d,d),1);
            MitX_h = arrayfun(@(ii) MitX_h(V(ii,:)), 1:n*m, 'un', 0);
            MitX_h = MitX_h';
            MitX_h = reshape(MitX_h,m,n);
            MitX_h = cell2mat(MitX_h);
            Mit = Mit + MitX_h;
        end
        
        MP = surf(X1,X2,Mit);
%         set(MP, 'edgecolor','none')
        set(MP,'LineStyle','none')
%         set(gca,'ZTickLabel',[])
    
    elseif (D==2.5) % 3D plot with conditional, mit is 1 dimensional
        x = -1:0.05:6;
        n = length(x);
        [X1,X2] = meshgrid(x,x);
        V1 = reshape(X1,n*n,1); V2 = reshape(X2,n*n,1);
        V = [V1,V2];
        Mit = zeros(n,m);

        for h=1:H
            p_h = mit.p(h);
            mu_h = mit.mu(h,:);
            Sigma_h = reshape(mit.Sigma(h,:),d,d);
            df_h = mit.df(h);

            MitX_h = @(a) p_h*dmvt(a, mu_h, Sigma_h, df_h, GamMat);
        %     MitX_h = @(a) 1*dmvt(a,mit_init.mu,reshape(mit_init.Sigma,d,d),1);
            MitX_h = arrayfun(@(ii) MitX_h(V(ii,:)), 1:n*m, 'un', 0);
            MitX_h = MitX_h';
            MitX_h = reshape(MitX_h,n,m);
            MitX_h = cell2mat(MitX_h);
            Mit = Mit + MitX_h';
        end
        
        for h = 1:H
            p_h = mit.p(h);
            mu_h = mit.mu(h,:);
            Sigma_h = mit.Sigma(h,:);
            df_h = mit.df(h);

            MitX_h = @(a) p_h*dmvt(a,mu_h,Sigma_h,df_h);
            MitX_h = arrayfun(@(ii) MitX_h(x(ii,:)), 1:n, 'un', 0);
            MitX_h = cell2mat(MitX_h);
            MitX_h = MitX_h';
            Mit = Mit + MitX_h;
        end
        
        surf(x,x,Mit);
%         set(gca,'ZTickLabel',[])
    end
end