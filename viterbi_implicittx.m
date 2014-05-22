function [path,totalcost,traceback] = viterbi_implicittx(lhoods, initialdist, relativetx)
% path = viterbi_implicittx(lhoods, initialdist, relativetx)
%   Find a viterbi path through a matrix of likelihood columns
%   where the transition matrix is described implicitly by a vector
%   of costs for relative bins.
%   relativetx is 2*nrows-1 long, and is shifted so that
%   relativetx(nrows) is always the "self-loop" prob for each
%   state. 
%   You want something like:
%     SM = simmx(DA, DB);
%     [nr,nc] = size(SM);
%     pri = exp(-0.5*([1:nr]/50).^2);
%     rtx = max(0.001, exp(-abs([-(nr-2):(nr)]'/2))); (or maybe 0.4
%                        instead of 2 for sharper diagonal constraints)
%     [pp,tc,trb] = viterbi_implicittx(SM, pri, rtx);
% 2014-03-13 Dan Ellis dpwe@ee.columbia.edu

[nr, nc] = size(lhoods);

totalcost = zeros(nr, nc);
traceback = zeros(nr, nc);

% Normalize the local match probs to max 1 in each frame
llhoods = log(lhoods./repmat(max(lhoods), size(lhoods,1),1));
% Normalize the tx probs to max 1
ltx = log(relativetx/max(relativetx));
%disp(['log tx(0)=',num2str(ltx(nr))]);

totalcost(:,1) = log(initialdist([1:nr])') + llhoods(:,1);

for i = 2:nc
  for j = 1:nr
%    if i > 900
%      disp([num2str([size(totalcost(:,i)), size(relativetx([nr-j+ ...
%                          [1:nr]]))])]);
%    end
    [cost, prev] = max(totalcost(:, i-1) ...
                       + ltx([nr-j+[1:nr]]));
    totalcost(j,i) = cost + llhoods(j,i);
    traceback(j,i) = prev;
  end
  if mod(i,100)==0
    disp(['col ', num2str(i)]);
  end
end

% traceback
path = zeros(1,nc);
[cost, path(nc)] = max(totalcost(:,nc));
disp(['best path log prob per step = ', num2str(cost/nc)]);
for i = nc:-1:2
  path(i-1) = traceback(path(i),i);
end
