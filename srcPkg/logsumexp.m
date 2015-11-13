function c = logsumexp(a, b)
% Computes log(exp(a) + exp(b)) for matrices a and b while guarding against
% overfow if exp(a) and or exp(b) are very large

d = max(a, b);
d(d == -Inf) = 0;

c = log(exp(a - d) + exp(b - d)) + d;