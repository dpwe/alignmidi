function [b,onsetenv,oesr,D,cumscore] = beat2(d,sr,startbpm,tightness,doplot)
% [b,onsetenv,oesr,D,cumscore] = beat(d,sr,startbpm,tightness,doplot)
%   b returns the times (in sec) of the beats in the waveform d, samplerate sr.
%   startbpm specifies the target tempo.  If it is a two-element
%   vector, it is taken as the mode of a tempo search window, with 
%   the second envelope being the spread (in octaves) of the
%   search, and the best tempo is calculated (with tempo.m).
%   tightness controls how tightly the start tempo is enforced
%   within the beat (default 6, larger = more rigid); if it is a 
%   two-element vector the second parameter is alpha, the strength 
%   of transition costs relative to local match (0..1, default 0.7).
%   doplot enables diagnostic plots; if it has two elements, they
%   are the time range (in sec) for the diagnostic plots.
%   onsetenv returns the raw onset detection envelope
%   D returns the mel-spectrogram, 
%   cumscore returns the per-frame cumulated dynamic-programming score.
% 2006-08-25 dpwe@ee.columbia.edu
% uses: localmax

%   Copyright (c) 2006 Columbia University.
% 
%   This file is part of LabROSA-coversongID
% 
%   LabROSA-coversongID is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License version 2 as
%   published by the Free Software Foundation.
% 
%   LabROSA-coversongID is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with LabROSA-coversongID; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
%   02110-1301 USA
% 
%   See the file "COPYING" for the text of the license.

if nargin < 3;   startbpm = 0; end
if nargin < 4;   tightness = 0; end
if nargin < 5;   doplot = 0; end

PLOTCUMSCORE = 0;  % leave 4th pane as tempo autoco

if length(startbpm) == 2
  temposd = startbpm(2);
  startbpm = startbpm(1);
else
  temposd = 0; 
end
if length(tightness) == 2
  alpha = tightness(2);
  tightness = tightness(1);
else
  alpha = 0;
end
if tightness == 0;  tightness = 400; end

% Have we been given an envelope (nonnegative waveform)
%if min(d) >= 0  % but onsetenv is HPF'd, so no longer nonneg
% just look for an unlikely audio SR, more likely envelope
if sr < 2000
  onsetenv = d;
  oesr = sr;
%  disp(['beat: treating input as onset strength envelope']);
else
  onsetenv = [];
end

% debug/plotting options
plotlims = [];
if length(doplot) > 1
  % specify zoom-in limits too
  plotlims = doplot;
  doplot = 1;
end
if doplot > 0;  debug = 1; else debug = 0; end

b = [];

% Select tempo search either with startbpm = 0 (means use defaults)
% or startbpm > 0 but temposd > 0 too (means search around startbpm)
% If onsetenv is empty, have to run tempo too to convert waveform
% to onsetenv, but we might not use the tempo it picks.
if startbpm == 0 | temposd > 0 | length(onsetenv) == 0

  if startbpm == 0
    tempomean = 240;
  else
    tempomean = startbpm;
  end

  if temposd == 0
    temposd = 1.0;
  end
  
  % Subfunction estimates global BPM; returns 'onset strength'
  % waveform onsetenv
  % If we were given an onsetenv as input, will use that
  [t,xcr,D,onsetenv,oesr,ff] = tempo2(d,sr,tempomean,temposd,debug);
  
  % tempo.m returns the top-2 BPM estimates; use faster one for
  % beat tracking
  usemax = 0;
  if (startbpm == 0 | temposd > 0)
    if usemax == 1
      startbpm = max(t([1 2]));
    else
      % try actual preferred tempo
      if t(3) > .5
        startbpm = t(1);
      else
        startbpm = t(2);
      end
    end
      
  end

  if debug == 1
    % plot the mel-specgram
    tt = [1:length(onsetenv)]/oesr;
    subplot(411)
    imagesc(tt,[1 40],D); axis xy
    ytk = get(gca,'YTick');
    set(gca,'YTickLabel',round(ff(ytk)));
    subplot(412)
    plot(tt,onsetenv);
  
    disp(['startbpm=',num2str(startbpm)]);
  end

end

% AGC on onsetenv
onsetenv = onsetenv/std(onsetenv);

% convert startbpm to startpd
startpd = round((60*oesr)/startbpm);
%disp(['startpd=',num2str(startpd)]);

pd = startpd;

% Smooth beat events
templt = exp(-0.5*(([-pd:pd]/(pd/32)).^2));
localscore = conv(templt,onsetenv);
%localscore = localscore(round(length(templt)/2)+[1:length(onsetenv)]);
localscore = localscore(pd+[1:length(onsetenv)]);
%imagesc(localscore)%%%%

% DP version:
% backlink(time) is index of best preceding time for this point
% cumscore(time) is total cumulated score to this point

backlink = zeros(1,length(localscore));
cumscore = zeros(1,length(localscore));

% search range for previous beat
prange = round(-2*pd):-round(pd/2);

% Skewed window
txwt = (-tightness*abs((log(prange/-pd)).^2));

starting = 1;
for i = 1:length(localscore)
  
  timerange = i + prange;
  
  % Are we reaching back before time zero?
  zpad = max(0, min(1-timerange(1),length(prange)));

  % Search over all possible predecessors and apply transition 
  % weighting
  scorecands = txwt + [zeros(1,zpad),cumscore(timerange(zpad+1:end))];
  % Find best predecessor beat
  [vv,xx] = max(scorecands);
  % Add on local score
  cumscore(i) = vv + localscore(i) - alpha;

  % special case to catch first onset
%  if starting == 1 & localscore(i) > 100*abs(vv)
  if starting == 1 & localscore(i) < 0.01*max(localscore);
    backlink(i) = -1;
  else
    backlink(i) = timerange(xx);
    % prevent it from resetting, even through a stretch of silence
    starting = 0;
  end
  
end

%%%% Backtrace

% Cumulated score is stabilized to lie in constant range, 
% so just look for one near the end that has a reasonable score
medscore = median(cumscore(localmax(cumscore)));
%maxscore = max(cumscore);
%bestendx = max(find(cumscore .* localmax(cumscore) > 0.75*maxscore));

bestendposs = find(cumscore .* localmax(cumscore) > 0.5*medscore);
bestendx = max(bestendposs);

b = bestendx;

while backlink(b(end)) > 0
  b = [b,backlink(b(end))];
end

b = fliplr(b);

%subplot(414); plot(b/oesr,localscore(b));

% use the smoothed version of the onset env
onsetenv = localscore;

% Actually choose start and end looking only on the beattimes
boe = localscore(b);
bwinlen = 5;
sboe = conv(hanning(bwinlen),boe);
sboe = sboe(floor(bwinlen/2)+1:length(boe));
thsboe = 0.5*sqrt(mean(sboe.^2));
% Keep only beats from first to last time that 
% smoothed beat onset times exceeds the threshold
b = b(min(find(sboe>thsboe)):max(find(sboe>thsboe)));

% return beat times in secs
b = b / oesr;

% Now done better above...
%% remove beats beyond last substantial beat 
%oethresh = 1.5*(mean(onsetenv.^2)^.5)
%b = b(b < (max(find(onsetenv > oethresh))+pd/2)/oesr);
%% .. and in the beginning
%b = b(b > (min(find(onsetenv < oethresh))-pd/2)/oesr);

% Debug visualization
if doplot == 1
  subplot(411)
  hold on;
  plot([b;b],[0;40]*ones(1,length(b)),'w');
  hold off;

  subplot(412)
  hold on;
  plot([b;b],[-2;5]*ones(1,length(b)),'g');
  hold off;
  ax = axis;
  ax([3 4]) = [-2 5];
  axis(ax);

  % redo 3rd pane as xcorr with templt
  subplot(413)
  tt = [1:length(localscore)]/oesr;
  plot(tt,localscore);
  hold on; plot([b;b],[min(localscore);max(localscore)]*ones(1,length(b)),'g'); hold off
  hold on; plot(tt(bestendposs),localscore(bestendposs),'or'); hold off
  ax = axis;
  ax([3 4]) = [-10 80];
  axis(ax);
   
  if PLOTCUMSCORE
     % 4th pane as cumscore
     subplot(414)
     tt = [1:length(localscore)]/oesr;
     ocumscore = cumscore - [0:length(cumscore)-1]*max(cumscore)/length(cumscore);
     plot(tt,ocumscore);
     hold on; plot([b;b],[min(ocumscore);max(ocumscore)]*ones(1,length(b)),'g'); hold off
  end
     
  if length(plotlims) > 0
    for i = 1:(3+PLOTCUMSCORE);
      subplot(4,1,i)
      ax = axis;
      ax([1 2]) = plotlims;
      axis(ax);
    end
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,xcr,D,onsetenv,oesr,ff] = tempo2(d,sr,tmean,tsd,debug)
% [t,xcr,D,onsetenv,oesr,ff] = tempo(d,sr,tmean,tsd,debug)
%    Estimate the overall tempo of a track for the MIREX McKinney
%    contest.  
%    d is the input audio at sampling rate sr.  tmean is the mode
%    for BPM weighting (in bpm) and tsd is its spread (in octaves).
%    onsetenv is an already-calculated onset envelope (so d is
%    ignored).  debug causes a debugging plot.
%    Output t(1) is the lower BPM estimate, t(2) is the faster,
%    t(3) is the relative weight for t(1) compared to t(2).
%    xcr is the windowed autocorrelation from which the BPM peaks were picked.
%    D is the mel-freq spectrogram
%    onsetenv is the "onset strength waveform", used for beat tracking
%    oesr is the sampling rate of onsetenv and D.
%    ff returns the center freqs of each row of D.
%
% 2006-08-25 dpwe@ee.columbia.edu
% uses: localmax, fft2melmx

%   Copyright (c) 2006 Columbia University.
% 
%   This file is part of LabROSA-coversongID
% 
%   LabROSA-coversongID is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License version 2 as
%   published by the Free Software Foundation.
% 
%   LabROSA-coversongID is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with LabROSA-coversongID; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
%   02110-1301 USA
% 
%   See the file "COPYING" for the text of the license.

if nargin < 3;   tmean = 110; end
if nargin < 4;   tsd = 0.9; end
if nargin < 5;   debug = 0; end

if sr < 2000
  % we were passed an onset env, not a waveform
  oesr = sr;
  onsetenv = d;
  %disp('data taken as onset envelope');
else
  onsetenv = [];
  
%  sro = 8000;
%  sro = 22050;
  % 2013-07-19: Calculate onset function using all the bandwidth
  % we're given
  sro = sr;
  % specgram: 256 bin @ 8kHz = 32 ms / 4 ms hop
  %swin = 256;
  swin = 2^round(log(0.016*sro)/log(2));
  %shop = 32;
  shop = round(0.004*sro);
  % mel channels
  nmel = 40;
  % sample rate for specgram frames (granularity for rest of processing)
  oesr = sro/shop;
end
  
% autoco out to 4 s
acmax = round(4*oesr);

D = 0;
ff = [];

if length(onsetenv) == 0
  % no onsetenv provided - have to calculate it

  % ensure mono
  if size(d,2) > 1
    d = mean(d,2);
  end
  
  % resample to target sampling rate?
  if (sr ~= sro)
    gg = gcd(sro,sr);
    d = resample(d,sro/gg,sr/gg);
    sr = sro;
  end

  D = specgram(d,swin,sr,swin,swin-shop);

  % Construct db-magnitude-mel-spectrogram
%  mlmx = fft2melmx(swin,sr,nmel);
  [mlmx,ff] = fft2melmx(swin,1.0*sr,nmel);
  D = 20*log10(max(1e-10,mlmx(:,1:(swin/2+1))*abs(D)));

  % Only look at the top 80 dB
  D = max(D, max(max(D))-80);

  %imgsc(D)
  
  % The raw onset decision waveform
  mm = (mean(max(0,diff(D')')));
  eelen = length(mm);

  % dc-removed mm
  onsetenv = filter([1 -1], [1 -.99],mm);

end  % of onsetenv calc block

% Find rough global period
% Only use the 1st 90 sec to estimate global pd (avoid glitches?)

maxd = 60;
maxt = 120; % sec
maxcol = min(round(maxt*oesr),length(onsetenv));
mincol = max(1,maxcol-round(maxd*oesr));

xcr = xcorr(onsetenv(mincol:maxcol),onsetenv(mincol:maxcol),acmax);

% find local max in the global ac
rawxcr = xcr(acmax+1+[0:acmax]);

% window it around default bpm
bpms = 60*oesr./([0:acmax]+0.1);
xcrwin = exp(-.5*((log(bpms/tmean)/log(2)/tsd).^2));

xcr = rawxcr.*xcrwin;

xpks = localmax(xcr);  
% will not include any peaks in first down slope (before goes below
% zero for the first time)
xpks(1:min(find(xcr<0))) = 0;
% largest local max away from zero
maxpk = max(xcr(xpks));

% ?? then period is shortest period with a peak that approaches the max
%maxpkthr = 0.4;
%startpd = -1 + min(find( (xpks.*xcr) > maxpkthr*maxpk ) );
%startpd = -1 + (find( (xpks.*xcr) > maxpkthr*maxpk ) );

%% no, just largest peak after windowing
%startpd = -1 + find((xpks.*xcr) == max(xpks.*xcr));
%% ??Choose acceptable peak closest to 120 bpm
%%[vv,spix] = min(abs(60./(startpd/oesr) - 120));
%%startpd = startpd(spix);
%% No, just choose shortest acceptable peak
%startpd = startpd(1);
%
%% Choose best peak out of .33 .5 2 3 x this period
%candpds = round([.33 .5 2 3]*startpd);
%candpds = candpds(candpds < acmax);
%
%[vv,xx] = max(xcr(1+candpds));
%
%startpd2 = candpds(xx);


%% Add in 2x, 3x, choose largest combined peak
%xcr2 = resample(xcr,1,2);
%xcr2 = xcr2 + xcr(1:length(xcr2));
%xcr3 = resample(xcr,1,3);
%xcr3 = xcr3 + xcr(1:length(xcr3));
% Quick and dirty explicit downsampling
lxcr = length(xcr);
xcr00 = [0, xcr, 0];
%wts = exp(-wt^2);
%sc = 1/(1+2*wts);
%xcr2 = xcr(1:ceil(lxcr/2))+sc*(wts*xcr00(1:2:lxcr)+xcr00(2:2:lxcr+1)+wts*xcr00(3:2:lxcr+2));
%xcr3 = xcr(1:ceil(lxcr/3))+sc*(wts*xcr00(1:3:lxcr)+xcr00(2:3:lxcr+1)+wts*xcr00(3:3:lxcr+2));
xcr2 = xcr(1:ceil(lxcr/2))+.5*(.5*xcr00(1:2:lxcr)+xcr00(2:2:lxcr+1)+.5*xcr00(3:2:lxcr+2));
xcr3 = xcr(1:ceil(lxcr/3))+.33*(xcr00(1:3:lxcr)+xcr00(2:3:lxcr+1)+xcr00(3:3:lxcr+2));

%subplot(413)
%plot(xcr2);
%hold on;
%plot(xcr3,'c');
%hold off

if max(xcr2) > max(xcr3)
  [vv, startpd] = max(xcr2);
  startpd = startpd -1;
  startpd2 = startpd*2;
else
  [vv, startpd] = max(xcr3);
  startpd = startpd -1;
  startpd2 = startpd*3;
end

% Weight by superfactors
pratio = xcr(1+startpd)/(xcr(1+startpd)+xcr(1+startpd2));

t = [60/(startpd/oesr) 60/(startpd2/oesr) pratio];

% ensure results are lowest-first
if t(2) < t(1)
  t([1 2]) = t([2 1]);
  t(3) = 1-t(3);
end  

startpd = (60/t(1))*oesr;
startpd2 = (60/t(2))*oesr;

%  figure
%  disp(['tmean=',num2str(tmean),' tsd=',num2str(tsd),' maxpk=',num2str(startpd)]);
%  subplot(211)
%  plot([0:acmax],xcrwin/max(abs(xcrwin)),[0:acmax],xcr/max(abs(xcr)),...
%       [startpd startpd],[-1 1],'-r',[startpd2 startpd2],[-1 1],'-c')
%  subplot(212)
%  bpms(1) = bpms(2);
%  plot(bpms,xcrwin/max(abs(xcrwin)),bpms,xcr/max(abs(xcr)),...
%       [t(1) t(1)],[-1 1],'-r',[t(2) t(2)],[-1 1],'-c')

if debug > 0

  % Report results and plot weighted autocorrelation with picked peaks
  disp(['Global bt pd = ',num2str(t(1)),' @ ',num2str(t(3)),[' / ' ...
                      ''],num2str(t(2)),' bpm @ ',num2str(1-t(3))]);

  subplot(414)
  tt = [0:acmax]/oesr;
  plot(tt,xcr,'-b', ...
       tt,xcrwin*maxpk,'-r', ...
       [startpd startpd]/oesr, [min(xcr) max(xcr)], '-g', ...
       [startpd2 startpd2]/oesr, [min(xcr) max(xcr)], '-c');
  grid;

end

% Read in all the tempo settings
% for i = 1:20; f = fopen(['mirex-beattrack/train/train',num2str(i),'-tempo.txt']); r(i,:) = fscanf(f, '%f\n'); fclose(f); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = localmax(x)
% return 1 where there are local maxima in x (columnwise).
% don't include first point, maybe last point

[nr,nc] = size(x);

if nr == 1
  lx = nc;
elseif nc == 1
  lx = nr;
  x = x';
else
  lx = nr;
end

if (nr == 1) || (nc == 1)

  m = (x > [x(1),x(1:(lx-1))]) & (x >= [x(2:lx),1+x(lx)]);

  if nc == 1
    % retranspose
    m = m';
  end
  
else
  % matrix
  lx = nr;
  m = (x > [x(1,:);x(1:(lx-1),:)]) & (x >= [x(2:lx,:);1+x(lx,:)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wts,binfrqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
% [wts,frqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
%      Generate a matrix of weights to combine FFT bins into Mel
%      bins.  nfft defines the source FFT size at sampling rate sr.
%      Optional nfilts specifies the number of output bands required 
%      (else one per bark), and width is the constant width of each 
%      band relative to standard Mel (default 1).
%      While wts has nfft columns, the second half are all zero. 
%      Hence, Mel spectrum is fft2melmx(nfft,sr)*abs(fft(xincols,nfft));
%      minfrq is the frequency (in Hz) of the lowest band edge;
%      default is 0, but 133.33 is a common standard (to skip LF).
%      maxfrq is frequency in Hz of upper edge; default sr/2.
%      You can exactly duplicate the mel matrix in Slaney's mfcc.m
%      as fft2melmx(512, 8000, 40, 1, 133.33, 6855.5, 0);
%      htkmel=1 means use HTK's version of the mel curve, not Slaney's.
%      constamp=1 means make integration windows peak at 1, not sum to 1.
%      frqs returns bin center frqs.
% 2004-09-05  dpwe@ee.columbia.edu  based on fft2barkmx

if nargin < 2;     sr = 8000;      end
if nargin < 3;     nfilts = 40;    end
if nargin < 4;     width = 1.0;    end
if nargin < 5;     minfrq = 0;     end  % default bottom edge at 0
if nargin < 6;     maxfrq = sr/2;  end  % default top edge at nyquist
if nargin < 7;     htkmel = 0;     end
if nargin < 8;     constamp = 0;   end


wts = zeros(nfilts, nfft);

% Center freqs of each FFT bin
fftfrqs = [0:(nfft/2)]/nfft*sr;

% 'Center freqs' of mel bands - uniformly spaced between limits
minmel = hz2mel(minfrq, htkmel);
maxmel = hz2mel(maxfrq, htkmel);
binfrqs = mel2hz(minmel+[0:(nfilts+1)]/(nfilts+1)*(maxmel-minmel), htkmel);

binbin = round(binfrqs/sr*(nfft-1));

for i = 1:nfilts
%  fs = mel2hz(i + [-1 0 1], htkmel);
  fs = binfrqs(i+[0 1 2]);
  % scale by width
  fs = fs(2)+width*(fs - fs(2));
  % lower and upper slopes for all bins
  loslope = (fftfrqs - fs(1))/(fs(2) - fs(1));
  hislope = (fs(3) - fftfrqs)/(fs(3) - fs(2));
  % .. then intersect them with each other and zero
%  wts(i,:) = 2/(fs(3)-fs(1))*max(0,min(loslope, hislope));
  wts(i,1+[0:(nfft/2)]) = max(0,min(loslope, hislope));

  % actual algo and weighting in feacalc (more or less)
%  wts(i,:) = 0;
%  ww = binbin(i+2)-binbin(i);
%  usl = binbin(i+1)-binbin(i);
%  wts(i,1+binbin(i)+[1:usl]) = 2/ww * [1:usl]/usl;
%  dsl = binbin(i+2)-binbin(i+1);
%  wts(i,1+binbin(i+1)+[1:(dsl-1)]) = 2/ww * [(dsl-1):-1:1]/dsl;
% need to disable weighting below if you use this one

end

if (constamp == 0)
  % Slaney-style mel is scaled to be approx constant E per channel
  wts = diag(2./(binfrqs(2+[1:nfilts])-binfrqs(1:nfilts)))*wts;
end

% Make sure 2nd half of FFT is zero
wts(:,(nfft/2+1):nfft) = 0;
% seems like a good idea to avoid aliasing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = mel2hz(z, htk)
%   f = mel2hz(z, htk)
%   Convert 'mel scale' frequencies into Hz
%   Optional htk = 1 means use the HTK formula
%   else use the formula from Slaney's mfcc.m
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
  htk = 0;
end

if htk == 1
  f = 700*(10.^(z/2595)-1);
else
  
  f_0 = 0; % 133.33333;
  f_sp = 200/3; % 66.66667;
  brkfrq = 1000;
  brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
  logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

  linpts = (z < brkpt);

  f = 0*z;

  % fill in parts separately
  f(linpts) = f_0 + f_sp*z(linpts);
  f(~linpts) = brkfrq*exp(log(logstep)*(z(~linpts)-brkpt));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = hz2mel(f,htk)
%  z = hz2mel(f,htk)
%  Convert frequencies f (in Hz) to mel 'scale'.
%  Optional htk = 1 uses the mel axis defined in the HTKBook
%  otherwise use Slaney's formula
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
  htk = 0;
end

if htk == 1
  z = 2595 * log10(1+f/700);
else
  % Mel fn to match Slaney's Auditory Toolbox mfcc.m

  f_0 = 0; % 133.33333;
  f_sp = 200/3; % 66.66667;
  brkfrq = 1000;
  brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
  logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

  linpts = (f < brkfrq);

  z = 0*f;

  % fill in parts separately
  z(linpts) = (f(linpts) - f_0)/f_sp;
  z(~linpts) = brkpt+(log(f(~linpts)/brkfrq))./log(logstep);

end
