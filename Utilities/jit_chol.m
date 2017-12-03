% @author: Maziar Raissi

function L = jit_chol(K)

jitter = (1e-16)*abs(mean(diag(K)));
num_tries = 0;
max_tries = 16;
N = size(K,1);
[L, p] = chol(K + jitter*eye(N), 'lower');

while p > 0 && num_tries < max_tries
    jitter = 10*jitter;
    num_tries = num_tries + 1;
    [L, p] = chol(K + jitter*eye(N), 'lower');
end

if p > 0
    disp(jitter);
    fprintf(1,'Covariance is ill-conditioned\n');
end

end