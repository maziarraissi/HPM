function [K] = kv1v0(x, xp, hyp, ubarp, vbarp, dt, i)

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

ubarp = repmat(ubarp',n_x,1);
vbarp = repmat(vbarp',n_x,1);

switch i


case 0

K=exp(1).^(logsigmav+(-1/2).*exp(1).^((-1).*logthetav).*(x+(-1).*xp).^2); ...
  


case 1 % logsigmau

K=0;


case 2 % logthetau

K=0;


case 3 % logsigmav

K=exp(1).^(logsigmav+(-1/2).*exp(1).^((-1).*logthetav).*(x+(-1).*xp).^2); ...
  


case 4 % logthetav

K=(1/2).*exp(1).^(logsigmav+(-1).*logthetav+(-1/2).*exp(1).^((-1).* ...
  logthetav).*(x+(-1).*xp).^2).*(x+(-1).*xp).^2;


case 5 % a1

K=0;


case 6 % a2

K=0;


otherwise
        
        K = zeros(n_x, n_xp);
end

if K == 0

    K = zeros(n_x, n_xp);

end

end
