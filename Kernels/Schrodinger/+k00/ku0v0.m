function [K] = ku0v0(x, xp, hyp, ubar, vbar, ubarp, vbarp, dt, i)

logsigmau = hyp(1);
logthetau = hyp(2);
logsigmav = hyp(3);
logthetav = hyp(4);

a1 = hyp(5);
a2 = hyp(6);

n_x = size(x,1);
n_xp = size(xp,1);

x = repmat(x,1,n_xp);
xp = repmat(xp',n_x,1);

ubar = repmat(ubar,1,n_xp);
vbar = repmat(vbar,1,n_xp);
ubarp = repmat(ubarp',n_x,1);
vbarp = repmat(vbarp',n_x,1);

switch i


case 0

K=dt.*(a2.*exp(1).^(logsigmav+(-1/2).*exp(1).^((-1).*logthetav).*(x+(-1).* ...
  xp).^2).*(ubar.^2+vbar.^2)+(-1).*a2.*exp(1).^(logsigmau+(-1/2).*exp(1) ...
  .^((-1).*logthetau).*(x+(-1).*xp).^2).*(ubarp.^2+vbarp.^2)+a1.*exp(1).^( ...
  logsigmau+(-2).*logthetau+(-1/2).*exp(1).^((-1).*logthetau).*(x+(-1).* ...
  xp).^2).*(exp(1).^logthetau+(-1).*(x+(-1).*xp).^2)+(-1).*a1.*exp(1).^( ...
  logsigmav+(-2).*logthetav+(-1/2).*exp(1).^((-1).*logthetav).*(x+(-1).* ...
  xp).^2).*(exp(1).^logthetav+(-1).*(x+(-1).*xp).^2));


case 1 % logsigmau

K=(-1).*dt.*exp(1).^(logsigmau+(-2).*logthetau+(-1/2).*exp(1).^((-1).* ...
  logthetau).*(x+(-1).*xp).^2).*(a2.*exp(1).^(2.*logthetau).*(ubarp.^2+ ...
  vbarp.^2)+a1.*((-1).*exp(1).^logthetau+(x+(-1).*xp).^2));


case 2 % logthetau

K=(-1/2).*dt.*exp(1).^(logsigmau+(-3).*logthetau+(-1/2).*exp(1).^((-1).* ...
  logthetau).*(x+(-1).*xp).^2).*(a1.*(2.*exp(1).^(2.*logthetau)+(-5).*exp( ...
  1).^logthetau.*(x+(-1).*xp).^2+(x+(-1).*xp).^4)+a2.*exp(1).^(2.* ...
  logthetau).*(ubarp.^2+vbarp.^2).*(x+(-1).*xp).^2);


case 3 % logsigmav

K=dt.*exp(1).^(logsigmav+(-2).*logthetav+(-1/2).*exp(1).^((-1).*logthetav) ...
  .*(x+(-1).*xp).^2).*(a2.*exp(1).^(2.*logthetav).*(ubar.^2+vbar.^2)+a1.*( ...
  (-1).*exp(1).^logthetav+(x+(-1).*xp).^2));


case 4 % logthetav

K=(1/2).*dt.*exp(1).^(logsigmav+(-3).*logthetav+(-1/2).*exp(1).^((-1).* ...
  logthetav).*(x+(-1).*xp).^2).*(a1.*(2.*exp(1).^(2.*logthetav)+(-5).*exp( ...
  1).^logthetav.*(x+(-1).*xp).^2+(x+(-1).*xp).^4)+a2.*exp(1).^(2.* ...
  logthetav).*(ubar.^2+vbar.^2).*(x+(-1).*xp).^2);


case 5 % a1

K=dt.*(exp(1).^(logsigmau+(-2).*logthetau+(-1/2).*exp(1).^((-1).* ...
  logthetau).*(x+(-1).*xp).^2).*(exp(1).^logthetau+(-1).*(x+(-1).*xp).^2)+ ...
  exp(1).^(logsigmav+(-2).*logthetav+(-1/2).*exp(1).^((-1).*logthetav).*( ...
  x+(-1).*xp).^2).*((-1).*exp(1).^logthetav+(x+(-1).*xp).^2));


case 6 % a2

K=dt.*(exp(1).^(logsigmav+(-1/2).*exp(1).^((-1).*logthetav).*(x+(-1).*xp) ...
  .^2).*(ubar.^2+vbar.^2)+(-1).*exp(1).^(logsigmau+(-1/2).*exp(1).^((-1).* ...
  logthetau).*(x+(-1).*xp).^2).*(ubarp.^2+vbarp.^2));


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
