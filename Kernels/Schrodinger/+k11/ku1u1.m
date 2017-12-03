function [K] = ku1u1(x, xp, hyp, i)

logsigmau = hyp(1);
logthetau = hyp(2);
logsigmav = hyp(3);
logthetav = hyp(4);

n_x = size(x,1);
n_xp = size(xp,1);

x = repmat(x,1,n_xp);
xp = repmat(xp',n_x,1);

switch i


case 0

K=exp(1).^(logsigmau+(-1/2).*exp(1).^((-1).*logthetau).*(x+(-1).*xp).^2); ...
  


case 1 % logsigmau

K=exp(1).^(logsigmau+(-1/2).*exp(1).^((-1).*logthetau).*(x+(-1).*xp).^2); ...
  


case 2 % logthetau

K=(1/2).*exp(1).^(logsigmau+(-1).*logthetau+(-1/2).*exp(1).^((-1).* ...
  logthetau).*(x+(-1).*xp).^2).*(x+(-1).*xp).^2;


case 3 % logsigmav

K=0;


case 4 % logthetav

K=0;


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
