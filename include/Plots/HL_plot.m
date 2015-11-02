figure(103)
aa = 0.95:0.01:1.05; n = length(aa);
ee = -4:0.01:-2; m =length(ee);
[AA,EE] = meshgrid(aa,ee);

Val = sqrt(AA).*EE;
Val =  100*(exp(Val/100) - 1);
MP = surf(AA,EE,Val);
hold on
surf (AA,EE,VaR_prelim*ones(m,n))
hold of

%%
figure(104)

kernel = @(x) posterior_debug_hl(x, y, a, b, mean(VaR_prelim), false);

xx = 0.85:0.01:1.15;  n = length(xx);
yy = -5:0.01:5;       m = length(yy);
%     xx = 0.95:0.01:1.05;
%     yy = -4:0.01:-2;

[XX,YY] = meshgrid(xx,yy);
XX1 = reshape(XX,n*m,1);
YY1 = reshape(YY,n*m,1);

Val1 = kernel([XX1,YY1]);
Val = reshape(Val1,m,n);
P = surf(XX,YY,Val);

%%
figure(105)

kernel = @(a) posterior_arch_hl(a, data, S, VaR_prelim, false);

xx = 0:0.01:0.5; n = length(xx);
yy = -5:0.01:5;  m = length(yy);

[XX,YY] = meshgrid(xx,yy);
XX1 = reshape(XX,n*m,1);
YY1 = reshape(YY,n*m,1);

Val1 = kernel([XX1,YY1]);
Val = reshape(Val1,m,n);
P = surf(XX,YY,Val);