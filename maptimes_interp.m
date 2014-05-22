function u = maptimes_interp(t,intime,outtime)
% u = maptimes(t,intime,outtime)
%   map the times in t according to the mapping that each point 
%   in intime corresponds to that value in outtime
% 2008-03-20 Dan Ellis dpwe@ee.columbia.edu

[tr, tc] = size(t);
t = t(:)';  % make into a row
nt = length(t);
nr = length(intime);

% Make sure both time ranges start at or before zero
pregap = max(intime(1), outtime(1));
intime = [intime(1) - pregap, intime];
outtime = [outtime(1) - pregap, outtime];

% Make sure there's a point beyond the end of both sequences
din = diff([intime, intime(end)+1]);
dout = diff([outtime, outtime(end)+1]);

% Decidedly faster than outer-product-array way
u = t(:);
for i = 1:nt
  ix = -1+min([find(intime > t(i)),length(outtime)]);
  % offset from that time
  dt = t(i) - intime(ix);
  % perform linear interpolation
  u(i) = outtime(ix) + (dt./din(ix))*dout(ix);
end
u = reshape(u,tr,tc);

