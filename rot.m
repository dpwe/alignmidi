function y = rot(x, n)
% Y = rot(X,N)                            Rotates the columns of X by N points.
%	i.e. the vector that used to be X(:,N) is now X(:,1), and the 
%	point that used to be X(:,1) is now X(:,size(X)-(N-1))
% dpwe 1994may23
sx = size(x);
s = sx(2);	% we're rotating the cols
rows = sx(1);
if(length(n) == 1)
  n = rem(n+s, s);
  y = [x(:,(n+1):s), x(:,1:n)];
elseif (length(n) == rows)
  y = zeros(rows, s);
  for r=1:rows
    nn = n(r);
    nn = rem(nn+s,s);
    y(r,:) = [x(r,(nn+1):s), x(r,1:nn)];
  end
else
  error('second argument must be scalar, or one value for each row')
end
