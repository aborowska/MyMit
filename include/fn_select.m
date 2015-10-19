function ind = fn_select(w, pc)
% get the indeces of draws with the largest IS weights
% w - vector of IS weights
% pc - fraction of highest weighted IS draws to use
% ind - N vector of indeces with highest IS weight draws
    N = size(w,1);
    Nl = round(N*pc);
	[~,ind] = sortrows(w);
    ind = ind((N-Nl+1):N,:);
%     w_a = [w,(1:N)'];
%     w_a = sortrows(w_a,1);
%     w_a = w_a((N-Nl+1):N,:);
%     ind = w_a(:,2);
%     ind = sort(ind);
%     w_hight = w(ind);
end
