function X = beatavg(Y,bts)
% X = beatavg(Y,bys)
%    Calculate average of columns of Y according to grid defined 
%    (real-valued) column indices in vector bts.
%    For folding spectrograms down into beat-sync features.
% 2006-09-26 dpwe@ee.columbia.edu

% beat-based segments
%bts = beattrack(d,sr);
%bttime = mean(diff(bts));
% map beats to specgram slices
ncols = size(Y,2);
coltimes = 0:(ncols-1);
% Ensure last value is end of array, and it only occurs once
bts = bts(bts < max(coltimes));
btse = [bts,max(coltimes)];
nbts = length(bts);
cols2beats = zeros(nbts, ncols);
for b = 1:nbts
  cols2beats(b,:) = ((coltimes >= btse(b)) & (coltimes < btse(b+1)))*1/(btse(b+1)-btse(b));
end

% The actual desired output
X = Y * cols2beats';
