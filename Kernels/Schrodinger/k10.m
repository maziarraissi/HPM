function K = k10(x, xp, hyp, Ubarp, dt, i)

ubarp = Ubarp(:,1);
vbarp = Ubarp(:,2);

Ku1u0 = k10.ku1u0(x, xp, hyp, ubarp, vbarp, dt, i);
Ku1v0 = k10.ku1v0(x, xp, hyp, ubarp, vbarp, dt, i);
Kv1u0 = k10.kv1u0(x, xp, hyp, ubarp, vbarp, dt, i);
Kv1v0 = k10.kv1v0(x, xp, hyp, ubarp, vbarp, dt, i);

K = [Ku1u0 Ku1v0;
     Kv1u0 Kv1v0];

end