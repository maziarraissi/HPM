function K = k00(x, xp, hyp, Ubar, Ubarp, dt, i)

ubar = Ubar(:,1);
vbar = Ubar(:,2);

ubarp = Ubarp(:,1);
vbarp = Ubarp(:,2);

Ku0u0 = k00.ku0u0(x, xp, hyp, ubar, vbar, ubarp, vbarp, dt, i);
Ku0v0 = k00.ku0v0(x, xp, hyp, ubar, vbar, ubarp, vbarp, dt, i);
Kv0u0 = Ku0v0';
Kv0v0 = k00.kv0v0(x, xp, hyp, ubar, vbar, ubarp, vbarp, dt, i);

K = [Ku0u0 Ku0v0;
     Kv0u0 Kv0v0];

end