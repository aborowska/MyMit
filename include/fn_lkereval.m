function lnk = fn_lkereval(kernel, theta)
% lnk - N vector of log-kernel evaluations at draws theta
    [N,~] = size(theta);
%     lnk = arrayfun(@(ii) kernel(theta(ii,:)), 1:N, 'un', 0);
%     lnk = cell2mat(lnk);
%     lnk = lnk';
    lnk = zeros(N,1);
    for ii = 1:N
        lnk(ii,1) = kernel(theta(ii,:));
    end
end