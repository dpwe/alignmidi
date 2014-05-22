function Y = normftrcols(X,W)
% Y = normftrcols(X,W)
%        X is a matrix of features, with each column representing a 
%        different sample (e.g. spectrogram-style).  Normalize each 
%        column by estimating mean and variance over a W-point sliding
%        window (def 51), and subtracting/dividing them out.
% 2003-03-16 dpwe@ee.columbia.edu

if nargin < 2
  W = 51;
end

h51 = hann(W);
sh51 = sum(h51);
W2 = floor(W/2);
nr = size(X,1);

Y = 0*X;

for c = 1:size(X,2)
  xx = X(:,c);
  mxx = conv(xx,h51)/sh51; 
  mxx=mxx(W2+[1:nr]); 
  vxx=sqrt(conv(h51,(xx-mxx).^2)/sh51); 
  vxx=vxx(W2+[1:nr]);
  % Don't divide by zero
  vxx = vxx + (vxx==0);
  Y(:,c) = (xx-mxx)./vxx; 
end

