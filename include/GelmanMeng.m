function r = GelmanMeng(x, A, B, C1, C2, L)
  x1 = x(:,1); x2 = x(:,2);
  r = -0.5*(A*(x1.^2).*(x2.^2) + x1.^2 + x2.^2 - 2*B.*x1.*x2 - 2*C1.*x1 - 2*C2.*x2);
  if (L==false)
    r = exp(r);
  end
end
