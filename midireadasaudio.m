function [d,sr] = midireadasaudio(MF,sr,domono)
% [d,sr] = readmidiasaudio(MF,sr,domono)
%    Read a midi file into an audio waveform.
%    Uses midi2wav to convert.
% 2013-07-17 Dan Ellis dpwe@ee.columbia.edu

af = midi2wav(MF);
[d,sr] = audioread(af, sr, domono);
delete(af);
