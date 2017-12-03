function K = k11(X, Xp, hyp, i)

x = X(:,1);
y = X(:,2);

xp = Xp(:,1);
yp = Xp(:,2);

Ku1u1 = k11.ku1u1(x, y, xp, yp, hyp, i);
Ku1v1 = k11.ku1v1(x, y, xp, yp, hyp, i);
Kv1u1 = k11.kv1u1(x, y, xp, yp, hyp, i);
Kv1v1 = k11.kv1v1(x, y, xp, yp, hyp, i);

K = [Ku1u1 Ku1v1;
     Kv1u1 Kv1v1];

end