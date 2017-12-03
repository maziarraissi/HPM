function f = Hypergeometric1F1(a,b,x)

% This function estimates the Kummer function with the specified tolerance
% the generalized hypergeometric series, noted below.  This solves Kummer's
% differential equation:
%
%       x*g''(x) + (b - x)*g'(x) - a*g(x) = 0

% Default tolerance is tol = 1e-10.  Feel free to change this as needed.

% Estimates the value by summing powers of the generalized hypergeometric
% series:
%
%       sum(n=0-->Inf)[(a)_n*x^n/{(b)_n*n!}
%
% until the specified tolerance is acheived.

term = x*a/b;
f = 1 + term;
n = 1;
an = a;
bn = b;
nmin = 150;
while(n < nmin)
  n = n + 1;
  an = an + 1;
  bn = bn + 1;
  term = x.*term*an/bn/n;
  f = f + term;
end

% VERSION INFORMATION
% v1 - Written to support only scalar inputs for x
% v2 - Changed to support column inputs for x by using the repmat
% command and using matrix multiplication to achieve the desired sum
%
% v3 - Credit goes to Ben Petschel for making this suggestion.
%    The previous method of creating vectors for multiplication to
%    produce the sum was replaced by a while loop that executes
%    until a certain tolerance is achieved.  My previous thinking
%    was avoiding a loop would produce a code that would execute
%    faster.  Ben pointed out this is not necessarily true, and not
%    true in this case.  Not only does the while loop used execute
%    faster for this calculation, but it is also more accurate.

end