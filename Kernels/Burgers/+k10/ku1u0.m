function [K] = ku1u0(x, xp, hyp, ubarp, dt, i)

logsigma = hyp(1);
logtheta = hyp(2);

a1 = hyp(3);
a2 = hyp(4);

n_x = size(x,1);
n_xp = size(xp,1);

x = repmat(x,1,n_xp);
xp = repmat(xp',n_x,1);

ubarp = repmat(ubarp',n_x,1);

switch i


case 0

K=exp(1).^(logsigma+(-2).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*(x+( ...
  -1).*xp).^2).*(exp(1).^logtheta.*(exp(1).^logtheta+a1.*dt.*ubarp.*(x+( ...
  -1).*xp))+a2.*dt.*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2));


case 1 % logsigma

K=exp(1).^(logsigma+(-2).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*(x+( ...
  -1).*xp).^2).*(exp(1).^logtheta.*(exp(1).^logtheta+a1.*dt.*ubarp.*(x+( ...
  -1).*xp))+a2.*dt.*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2));


case 2 % logtheta

K=(1/2).*exp(1).^(logsigma+(-3).*logtheta+(-1/2).*exp(1).^((-1).*logtheta) ...
  .*(x+(-1).*xp).^2).*((-1).*a2.*dt.*(2.*exp(1).^(2.*logtheta)+(-5).*exp( ...
  1).^logtheta.*(x+(-1).*xp).^2+(x+(-1).*xp).^4)+exp(1).^logtheta.*((-1).* ...
  a1.*dt.*ubarp.*(2.*exp(1).^logtheta+(-1).*(x+(-1).*xp).^2)+exp(1) ...
  .^logtheta.*(x+(-1).*xp)).*(x+(-1).*xp));


case 3 % a1

K=dt.*exp(1).^(logsigma+(-1).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*( ...
  x+(-1).*xp).^2).*ubarp.*(x+(-1).*xp);


case 4 % a2

K=dt.*exp(1).^(logsigma+(-2).*logtheta+(-1/2).*exp(1).^((-1).*logtheta).*( ...
  x+(-1).*xp).^2).*(exp(1).^logtheta+(-1).*(x+(-1).*xp).^2);


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
