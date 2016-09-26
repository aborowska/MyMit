function cont =  MitISEM_Control
    cont.mit.N = 10000; %1e5;
    cont.mit.Hmax = 10;
%  
%     cont.mit.Hmax1 = 10;  
%     cont.mit.Hmax2 = 1;      
    
    cont.mit.CV_tol = 0.1; %0.1
    cont.mit.CV_old = 100;

    cont.mit.ISpc = 0.1;
    cont.mit.pnc = 0.1; % probability of a new component
    cont.mit.dfnc = 1; % degrees of freedom of a new component
    cont.mit.tol_pr = 0;

    cont.mit.norm = true;

    cont.EM.maxit = 100;
    cont.EM.tol = 0.001;
    cont.EM.stab = 1e5;

    cont.df.maxit = 1000;
    cont.df.opt = true;
    cont.df.range = [1,10];
    cont.df.tol = eps^(0.25);

    cont.resmpl_on = false;
    cont.mit.iter_max = 2;
end