function [K] = ku1u0(x, xp, hyp, ubarp, dt, i)

logsigma = hyp(1);
logtheta = hyp(2);

a1 = hyp(3);
a2 = hyp(4);
a3 = hyp(5);
a4 = hyp(6);

n_x = size(x,1);
n_xp = size(xp,1);

x = repmat(x,1,n_xp);
xp = repmat(xp',n_x,1);

ubarp = repmat(ubarp',n_x,1);

switch i


case 0

K=exp(1).^(logsigma+(-4).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*(x+( ...
  -1).*xp).^2).*(exp(1).^logtheta.*(exp(1).^(3.*logtheta)+(-1).*a2.*dt.* ...
  exp(1).^logtheta.*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2)+(-3).*a3.* ...
  dt.*exp(1).^logtheta.*(x+(-1).*xp)+a1.*dt.*exp(1).^(2.*logtheta).* ...
  ubarp.*(x+(-1).*xp)+a3.*dt.*(x+(-1).*xp).^3)+a4.*dt.*(3.*exp(1).^(2.* ...
  logtheta)+(-6).*exp(1).^logtheta.*(x+(-1).*xp).^2+(x+(-1).*xp).^4));


case 1 % logsigma

K=exp(1).^(logsigma+(-4).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*(x+( ...
  -1).*xp).^2).*(exp(1).^logtheta.*(exp(1).^(3.*logtheta)+(-1).*a2.*dt.* ...
  exp(1).^logtheta.*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2)+(-3).*a3.* ...
  dt.*exp(1).^logtheta.*(x+(-1).*xp)+a1.*dt.*exp(1).^(2.*logtheta).* ...
  ubarp.*(x+(-1).*xp)+a3.*dt.*(x+(-1).*xp).^3)+a4.*dt.*(3.*exp(1).^(2.* ...
  logtheta)+(-6).*exp(1).^logtheta.*(x+(-1).*xp).^2+(x+(-1).*xp).^4));


case 2 % logtheta

K=(1/2).*exp(1).^(logsigma+(-5).*logtheta+(-1/2).*exp(1).^((-1).*logtheta) ...
  .*(x+(-1).*xp).^2).*((-1).*a4.*dt.*(12.*exp(1).^(3.*logtheta)+(-39).* ...
  exp(1).^(2.*logtheta).*(x+(-1).*xp).^2+14.*exp(1).^logtheta.*(x+(-1).* ...
  xp).^4+(-1).*(x+(-1).*xp).^6)+exp(1).^logtheta.*(a2.*dt.*exp(1) ...
  .^logtheta.*(2.*exp(1).^(2.*logtheta)+(-5).*exp(1).^logtheta.*(x+(-1).* ...
  xp).^2+(x+(-1).*xp).^4)+(x+(-1).*xp).*(a3.*dt.*(12.*exp(1).^(2.* ...
  logtheta)+(-9).*exp(1).^logtheta.*(x+(-1).*xp).^2+(x+(-1).*xp).^4)+(-1) ...
  .*exp(1).^(2.*logtheta).*(a1.*dt.*ubarp.*(2.*exp(1).^logtheta+(-1).*(x+( ...
  -1).*xp).^2)+exp(1).^logtheta.*((-1).*x+xp)))));


case 3 % a1

K=dt.*exp(1).^(logsigma+(-1).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*( ...
  x+(-1).*xp).^2).*ubarp.*(x+(-1).*xp);


case 4 % a2

K=(-1).*dt.*exp(1).^(logsigma+(-2).*logtheta+(-1/2).*exp(1).^((-1).* ...
  logtheta).*(x+(-1).*xp).^2).*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2);


case 5 % a3

K=dt.*exp(1).^(logsigma+(-3).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*( ...
  x+(-1).*xp).^2).*((-3).*exp(1).^logtheta+(x+(-1).*xp).^2).*(x+(-1).*xp); ...
  


case 6 % a4

K=dt.*exp(1).^(logsigma+(-4).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*( ...
  x+(-1).*xp).^2).*(3.*exp(1).^(2.*logtheta)+(-6).*exp(1).^logtheta.*(x+( ...
  -1).*xp).^2+(x+(-1).*xp).^4);


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
