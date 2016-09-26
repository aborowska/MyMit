function [pmit, CV_mix, CV, iter, pmit_step2, pmit_step3, pmit_adapt] = PMitISEM2(draw0, lnk0, w0, kernel, fn_const_X, fn_input_X, partition, d, cont, GamMat)
%     tic
profile on
    SS = length(partition);
    N = cont.mit.N;
    Hmax = cont.mit.Hmax;      
    iter_max = cont.mit.iter_max;

    CV_tol = cont.mit.CV_tol;
    CV_old = cont.mit.CV_old;
    CV = cell(SS,1);
    
 %% Step 0: Initialization
    % get/define initial mit density shape scale and degrees of freedom
    % aka: naive proposal density
    if isa(draw0, 'function_handle')
        mu_init = lnk0;
        kernel_init = draw0;
        k = size(mu_init,2);
        
        [mu, Sigma] = fn_initopt(kernel_init, mu_init);

        try
            r = fn_testSigma(Sigma);
        catch
            r = 1;
        end

        if (r == 1)
            Sigma = ones(k,k) + 2*diag(ones(k,1));
            Sigma = reshape(Sigma,1,k^2);
            display('Initial optimzation FAILED.')        
        else
            display('Initial optimzation OK.')
        end
        mit_init.mu = mu;
        mit_init.Sigma = Sigma;
        mit_init.df = cont.mit.dfnc;
        mit_init.p = 1;
    
        % get draws and IS weights from naive  
        [draw0, lnk0, ~] = fn_rmvgt_robust(N, mit_init, kernel, false);
        lnd0 = dmvgt(draw0, mit_init, true, GamMat);
        w0 = fn_ISwgts(lnk0, lnd0, cont.mit.norm);
    end
    
  CV0 =  fn_CVstop(w0, [], []);   
  tic
    %% STEP 1: ADAPTATION
    % Adapted patial mixture - to have an input to ISEM
    pmit_adapt = struct('mu',cell(1,SS),'Sigma',cell(1,SS),'df',cell(1,SS),'p',cell(1,SS));
    
    input_X = fn_input_X(draw0);
    for s = 1:SS    
        fprintf('\nStep 1. Subset no.: %d.\n',s)
        [s1, s2] = fn_partition_ends(partition, d, s);
        theta_s = draw0(:,s1:s2);
        if (s==1)
            [pmit_adapt(s).mu, pmit_adapt(s).Sigma] = fn_muSigma(theta_s, w0);
        else
            if ~isstruct(input_X)
                input_X = draw0(:,1:s1-1);
            else
                input_X.theta = draw0(:,1:s1-1);
            end              
            X = fn_const_X(input_X);
            [beta, Sigma] = fn_beta(theta_s,w0,X); 
            pmit_adapt(s).mu = beta; % x = theta(:,s1:s2);
            pmit_adapt(s).Sigma = Sigma;
        end   
        pmit_adapt(s).df = cont.mit.dfnc;
        pmit_adapt(s).p = 1;
    end


    if ~isstruct(input_X)
        input_X = draw0;
    else
        input_X.theta = draw0;
    end   
    lnk_adapt = kernel(draw0); 
    lnd_adapt = fn_dpmit2(input_X, pmit_adapt, partition, fn_const_X, true, GamMat);        
    w_adapt = fn_ISwgts(lnk_adapt, lnd_adapt, false); 
    CV_adapt = fn_CVstop(w_adapt, [], []);
      
    CV_mix = [CV0, CV_adapt];
    CV(:) = {CV0};
    
    pmit = pmit_adapt;
    hstop_mix = false;
    iter = 0;

%     while ((iter < iter_max) && (hstop_mix == false))
        iter = iter + 1;
        
        pmit_old = pmit;
        
        fprintf('\nPMit mixture iteration no.: %d.\n',iter)
        
        %% STEP 2: ITERATE ON NUMBER OF COMPONENTS
        for s = 1:SS
            fprintf('\nStep 2. Subset no.: %d.\n',s)

            hstop = false;
            H_s = length(pmit(s).p);  

            % the first component
            [s1, s2] = fn_partition_ends(partition, d, s);
            theta_s = draw0(:,s1:s2);
            if (s==1)
                pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat);
            else
                if ~isstruct(input_X)
                    input_X = draw0(:,1:s1-1);
                else
                    input_X.theta = draw0(:,1:s1-1);
                end              
                X = fn_const_X(input_X);
                if (iter == 1)
                    [beta, Sigma] = fn_beta(theta_s,w0,X); 
                    pmit(s).mu = beta; % x = theta(:,s1:s2);
                    pmit(s).Sigma = Sigma;
                end
                pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
            end
            
            if ~isstruct(input_X)
                input_X = draw0;
            else
                input_X.theta = draw0;
            end 
            lnd_curr = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat);

            w_curr = fn_ISwgts(lnk0, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
            CV_old = CV{s};
            [CV_new, ~] = fn_CVstop(w_curr, [], []);
            CV(s) = {[CV{s}, CV_new]}; 
            
            % Step 2b & 2c:
            while ((H_s < Hmax) && (hstop == false))
                H_s = H_s+1;
                ind_w = fn_select(w_curr, cont.mit.ISpc);
                theta_nc = draw0(ind_w,:);
%                 if isstruct(input_X)
%                     input_X_nc = structfun(@(xx) xx(ind_w,:), input_X, 'UniformOutput', false);
%                 end                
                w_nc = w_curr(ind_w);
                theta_nc_s = theta_nc(:,s1:s2);

                % NEW COMPONENT
                % compute new component's mode and scale from IS weights
                mit_nc.p = cont.mit.pnc;
                mit_nc.df = cont.mit.dfnc;
                if (s==1)
                    [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc_s, w_nc);
                else
%                     if ~isstruct(input_X)
%                         input_X_nc = theta_nc(:,1:s1-1);
%                     else
%                         input_X_nc.theta = theta_nc(:,1:s1-1);
%                     end 
                    X_nc = X(ind_w,:);
                    [beta, Sigma] = fn_beta(theta_nc_s,w_nc,X_nc); 
%                     X = fn_const_X(input_X_nc);
%                     X = fn_const_X(theta_nc(:,1:s1-1));
%                     [beta, Sigma] = fn_beta(theta_nc_s,w_nc,X); 
                    mit_nc.mu = beta; 
                    mit_nc.Sigma = Sigma;
                end              

                % COBINE OLD AND NC
                % combine the old mixture mit_new and the new component mit_nc
                pmit_old_s = pmit(s);
                pmit(s) = fn_updateMit(pmit(s), mit_nc); 

                % UPDATE COMBINED
                % update mode, scale and df  of all mixture components
                if (s==1)
                    pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat); % w_curr only to locate the new component; w0 to run ISEM
                else
%                     if ~isstruct(input_X)
%                         input_X = draw0(:,1:s1-1);
%                     else
%                         input_X.theta = draw0(:,1:s1-1);
%                     end              
%                     X = fn_const_X(input_X);
                    pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
                end
                
                if ~isstruct(input_X)
                    input_X = draw0;
                else
                    input_X.theta = draw0;
                end 
                lnd_curr = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat); % for all draws!

                w_curr = fn_ISwgts(lnk0, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
                CV_old = CV{s}; 
                CV_old = CV_old(:,end);
                [CV_new, hstop_new] = fn_CVstop(w_curr, CV_old, CV_tol);
                CV(s) = {[CV{s}, CV_new]}; 
        %         if (H > 1)
                    hstop = hstop_new;
        %         end   
            end
        end
        pmit_step2 = pmit;
time_pmit(1,1) = toc         
tic
%         Step 2 up
        for s = 1:SS 
            fprintf('\nStep 2 up. Subset no.: %d.\n',s)
            [s1, s2] = fn_partition_ends(partition, d, s);
            theta_s = draw0(:,s1:s2);
            if (s==1)
                pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat);
            else
                if ~isstruct(input_X)
                    input_X = draw0(:,1:s1-1);
                else
                    input_X.theta = draw0(:,1:s1-1);
                end              
                X = fn_const_X(input_X);
                pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
            end
        end
        pmit_step2_up = pmit;
time_step2_up =  toc 

% tic
        %% STEP 3: sample from pmit and check convergence
        [draw_pmit, lnk_pmit] = fn_p_rmvgt2(size(draw0,1), pmit, d, partition, kernel, fn_const_X, fn_input_X);         
        input_X = fn_input_X(draw_pmit);
%         if ~isstruct(input_X)
%             input_X = draw_pmit;
%         else
%             input_X.theta = draw_pmit;
%         end 
        lnd_pmit = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat);        
        
        w_pmit = fn_ISwgts(lnk_pmit, lnd_pmit, false); 
        CV_pmit = fn_CVstop(w_pmit, [], []);
        CV_mix = [CV_mix; CV_pmit, 0];
        fprintf('CV_mix(iter) = %6.4f ',CV_mix(end,:)); fprintf('\n');
        
        % update with the latest draws and the corresponding IS weighs
        for s = 1:SS 
            fprintf('\nStep 3. Subset no.: %d.\n',s)
            [s1, s2] = fn_partition_ends(partition, d, s);
            theta_s = draw_pmit(:,s1:s2);
            if (s==1)
                pmit(s) = fn_optimt(theta_s, pmit(s), w_pmit, cont, GamMat);
            else
                if ~isstruct(input_X)
                    input_X = draw_pmit(:,1:s1-1);
                else
                    input_X.theta = draw_pmit(:,1:s1-1);
                end              
                X = fn_const_X(input_X);
                pmit(s) = fn_Poptimt(theta_s, pmit(s), w_pmit, cont, GamMat, X);
            end
        end
        pmit_step3 = pmit;

        if ~isstruct(input_X)
            input_X = draw_pmit;
        else
            input_X.theta = draw_pmit;
        end 
        lnd_curr = fn_dpmit2(input_X, pmit, partition, fn_const_X, true, GamMat);  % for all draws!
        w_curr = fn_ISwgts(lnk_pmit, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
        [CV_curr, hstop_mix] = fn_CVstop(w_curr, CV_mix(end-1,2), CV_tol);
        CV_mix(end,2) = CV_curr;
        fprintf('CV_mix(iter) = %6.4f ',CV_mix(end,:)); fprintf('\n');

        draw0 = draw_pmit;
        w0 = w_pmit;
        lnk0 = lnk_pmit;
%  time_step3 = toc; 
%     end
%     if (CV_mix(end,2) > CV_mix(end-1,2))
%         pmit = pmit_old;
%     end
profile off
profile viewer
end