function par_KFS = IS_Model(par_NAIS, par_SV)
% Set the parameters for the artificial linear Gaussian model
% To be used in KFS recursions
    
    C = par_NAIS.C;
%     [n,S] = size(C);
%     Un = ones(n,S);
    
    c = par_SV(:,1);
    phi = par_SV(:,2);
    sigma2 = par_SV(:,3);
    P1 = sigma2./(1-phi.^2);
 
    par_KFS.P1 = P1';           % <-- VECTOR
    par_KFS.c = c';             % <-- VECTOR
%     par_KFS.c = c.*Un;
    par_KFS.H = C.^(-1);        % <-- MATRIX
%     par_KFS.Q = sigma2*Un;
%     par_KFS.d = 0*Un;
    par_KFS.Q = sigma2';        % <-- VECTOR
    par_KFS.d = 0;
%     par_KFS.T = phi*Un;
%     par_KFS.R = Un;
%     par_KFS.Z = Un;
    par_KFS.T = phi';           % <-- VECTOR
    par_KFS.R = 1;
    par_KFS.Z = 1;
end
