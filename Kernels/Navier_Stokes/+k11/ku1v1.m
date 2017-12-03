function [K] = ku1v1(x, y, xp, yp, hyp, i)

logsigma = hyp(1);
logthetax = hyp(2);
logthetay = hyp(3);

n_x = size(x,1);
n_y = size(y,1);
n_xp = size(xp,1);
n_yp = size(yp,1);

x = repmat(x,1,n_xp);
y = repmat(y,1,n_yp);
xp = repmat(xp',n_x,1);
yp = repmat(yp',n_y,1);

switch i


case 0

K=exp(1).^(logsigma+(-1).*logthetax+(-1).*logthetay+(-1/2).*exp(1).^((-1) ...
  .*logthetax).*(x+(-1).*xp).^2+(-1/2).*exp(1).^((-1).*logthetay).*(y+(-1) ...
  .*yp).^2).*(x+(-1).*xp).*(y+(-1).*yp);


case 1 % logsigma

K=exp(1).^(logsigma+(-1).*logthetax+(-1).*logthetay+(-1/2).*exp(1).^((-1) ...
  .*logthetax).*(x+(-1).*xp).^2+(-1/2).*exp(1).^((-1).*logthetay).*(y+(-1) ...
  .*yp).^2).*(x+(-1).*xp).*(y+(-1).*yp);


case 2 % logthetax

K=exp(1).^(logsigma+(-1).*logthetax+(-1).*logthetay+(-1/2).*exp(1).^((-1) ...
  .*logthetax).*(x+(-1).*xp).^2+(-1/2).*exp(1).^((-1).*logthetay).*(y+(-1) ...
  .*yp).^2).*((-1)+(1/2).*exp(1).^((-1).*logthetax).*(x+(-1).*xp).^2).*(x+ ...
  (-1).*xp).*(y+(-1).*yp);


case 3 % logthetay

K=exp(1).^(logsigma+(-1).*logthetax+(-1).*logthetay+(-1/2).*exp(1).^((-1) ...
  .*logthetax).*(x+(-1).*xp).^2+(-1/2).*exp(1).^((-1).*logthetay).*(y+(-1) ...
  .*yp).^2).*(x+(-1).*xp).*((-1)+(1/2).*exp(1).^((-1).*logthetay).*(y+(-1) ...
  .*yp).^2).*(y+(-1).*yp);


otherwise
        
        K = zeros(n_x, n_xp);
end

end
