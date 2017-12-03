function K = k00(X, Xp, hyp, Ubar, Ubarp, dt, i)

x = X(:,1);
y = X(:,2);

xp = Xp(:,1);
yp = Xp(:,2);

ubar = Ubar(:,1);
vbar = Ubar(:,2);

ubarp = Ubarp(:,1);
vbarp = Ubarp(:,2);

Ku0u0 = k00.ku0u0(x, y, xp, yp, hyp, ubar, vbar, ubarp, vbarp, dt, i);
Ku0v0 = k00.ku0v0(x, y, xp, yp, hyp, ubar, vbar, ubarp, vbarp, dt, i);
Kv0u0 = Ku0v0';
Kv0v0 = k00.kv0v0(x, y, xp, yp, hyp, ubar, vbar, ubarp, vbarp, dt, i);

K = [Ku0u0 Ku0v0;
     Kv0u0 Kv0v0];

end