function [d,sr] = midireadasaudio(MF,TARGETSR,FORCEMONO,START,DUR)
% [d,sr] = midireadasaudio(MF,TARGETSR,FORCEMONO,START,DUR)
%    Read a midi file into an audio waveform.
%    Uses midi2wav to convert.
%    Other args as audioread.
% 2013-07-17 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; TARGETSR = []; end
if nargin < 3; FORCEMONO = 0; end
if nargin < 4; START = 0; end
if nargin < 5; DUR = 0; end

af = midi2wav(MF);
[d,sr] = audioread(af, TARGETSR, FORCEMONO, START, DUR);
delete(af);
