function [p,q,D,phi,score] = dpmod(M,pen,G)
% [p,q,C,phi,score] = dp(M,pen,G) 
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q.
%    C returns total cost matrix.
%    pen is cost scale for for (0,1) and (1,0) steps (default 1.0)
%    G defines acceptable "gullies" (0..1, default 0).
%
% 2003-03-15 dpwe@ee.columbia.edu
% $Header: /Users/dpwe/projects/dtw/RCS/dp.m,v 1.2 2006/01/18 20:07:23 dpwe Exp $
% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

% penalty for horiz/vert 
if nargin < 2
  pen = 1;
end
if nargin < 3
  G = 0;
end

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
% add "gulleys" at start too
D(1,1:round(G*c)) = 0;
D(1:round(G*r),1) = 0;
% copy in local costs
D(2:(r+1), 2:(c+1)) = M;

% traceback
phi = zeros(r,c);

for i = 1:r;
  for j = 1:c;
    [dmin, tb] = min([D(i, j), pen+D(i, j+1), pen+D(i+1, j)]);
    D(i+1,j+1) = D(i+1,j+1)+dmin;
    phi(i,j) = tb;
  end
end

% Strip off the edges of the D matrix before returning
D = D(2:(r+1),2:(c+1));

if G == 0
  % Traceback from top left
  i = r; 
  j = c;
else
  % Traceback from lowest cost "to edge" (gulleys)
  TE = D(r,:);
  RE = D(:,c);
  % GulleyFix1: normalize points in gulleys by their shorter
  % theoretical lengths
  TE = TE ./ [1:length(TE)];
  RE = RE ./ [1:length(RE)]';
  % eliminate points not in gulleys
  TE(1:round((1-G)*c)) = NaN; %max(max(D));
  RE(1:round((1-G)*r)) = NaN; %max(max(D));
  if (min(TE) < min(RE))
    i = r;
    j = max(find(TE==min(TE)));
  else
    i = max(find(RE==min(RE)));
    j = c;
  end
end

score = D(i,j);
%disp(['dpmod: Best cost = ',num2str(score)]);

p = i;
q = j;

while i > 1 && j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else
    error;
  end
  p = [i,p];
  q = [j,q];
end

DEBUG = 0;
if DEBUG
  figure(2)
  subplot(321)
  plot([TE, fliplr(RE')])
  endpos = q(end);
  if endpos == c
    endpos = length(TE) + length(RE) - (p(end)-1);
  end
  ax = axis(); hold on; plot([endpos endpos], [ax(3) ax(4)], '-r'); ...
       hold off
  subplot(325)
  plot(p, q-p)
  subplot(122)
  imgsc(M)
  hold on; plot(q,p,'r'); hold off
  figure(1)
end
