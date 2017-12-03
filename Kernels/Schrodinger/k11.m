function K = k11(x, xp, hyp, i)

n_x = size(x,1);
n_xp = size(xp,1);

Ku1u1 = k11.ku1u1(x, xp, hyp, i);
Ku1v1 = zeros(n_x,n_xp);
Kv1u1 = zeros(n_x,n_xp);
Kv1v1 = k11.kv1v1(x, xp, hyp, i);

K = [Ku1u1 Ku1v1;
     Kv1u1 Kv1v1];

end