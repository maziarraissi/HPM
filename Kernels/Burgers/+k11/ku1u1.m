function [K] = ku1u1(x, xp, hyp, i)

logsigma = hyp(1);
logtheta = hyp(2);

n_x = size(x,1);
n_xp = size(xp,1);

x = repmat(x,1,n_xp);
xp = repmat(xp',n_x,1);

switch i


case 0

K=exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta).*(x+(-1).*xp).^2);


case 1 % logsigma

K=exp(1).^(logsigma+(-1/2).*exp(1).^((-1).*logtheta).*(x+(-1).*xp).^2);


case 2 % logtheta

K=(1/2).*exp(1).^(logsigma+(-1).*logtheta+(-1/2).*exp(1).^((-1).*logtheta) ...
  .*(x+(-1).*xp).^2).*(x+(-1).*xp).^2;


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
