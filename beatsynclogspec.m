function [D,b,frqs] = beatsynclogspec(d, sr, initbpm, tightness)
% [D,b,frqs] = beatsynclogspec(d, sr, initbpm, tightness)
%   Calculate beat-synchronous log-spectrogram features
%   D returns energy, b returns times of beats, 
%   frqs is center frq in Hz of each bin
% 2013-07-17 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; initbpm = [240 0.5]; end  % 2014-02-16 added .75 spread
if nargin < 4; tightness = 800; end

b = beat2(d, sr, initbpm, tightness);

% Make logspec features, at lower SR
d = resample(d, 1, 2);
sr = sr/2;

whitenord = 8;
twin = 0.050;
thop = 0.010;

fftlen = 2^round(log(twin*sr)/log(2));
hop = round(thop*sr);

fmin = 41.2; % Hz of low E on bass guitar
bpo = 12; % bins per octave

[Df,mx,frqs] = logfsgram(whiten(d,whitenord), ...
                         fftlen, sr, fftlen, fftlen-hop, fmin, bpo);

% Beat-sampled features
D = beatavg(Df, b/thop);
