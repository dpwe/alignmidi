%% alignmidi - Aligning MIDI files to music audio
%
% To obtain ground-truth transcriptions of real music audio, it is 
% sometimes possible to find a MIDI version of the track which can 
% then be aligned in time.  The time-aligned MIDI events can then 
% be taken as an approximate transcription of the audio.  There are 
% several ways to do this, including several choices of the domain 
% in which to do the matching; one popular approach is to
% synthesize the MIDI version to audio, then do a time alignment
% between the original audio and the synthetic audio (see, e.g.,
% Turetsky and Ellis ISMIR 2003, or Hu, Dannenberg, and Tzanetakis,
% WASPAA, 2003). 
%
% This matlab code refines that approach by first beat tracking both the 
% original audio and the resynthesized MIDI, then doing the time 
% warp alignment over units of complete beats.  Since the remapping 
% of MIDI event times will be done using each beat as an anchor
% point, the temporal precision is usually good (as long as the 
% beat tracker worked.
%
% Note that this code has been significantly improved in March
% 2014.  Firstly, we made a range of improvements to the local
% similarity computation and DTW path picking (v0.03).  Then we
% replaced DTW with full Viterbi alignment, which allows arbitrary
% jumps to accommodate missing or repeated sections etc.  You can
% still select DTW paths (which avoid really big jumps) as an
% option to alignmidi (params.usedtw=1); DTW runs a lot faster.

%% Example usage
% The code below shows how to run the alignment using the main
% <alignmidi.m alignmidi> routine, and a couple of ways of
% how the outputs can be used. 

% We will use a modified version of Christine Smit's readmidi_java
% from http://www.ee.columbia.edu/~csmit/karaoke_midi.html .
% If this doesn't work, you need to make sure KaraokeMidJava.jar
% is in Matlab's classpath.txt
if exist('readmidi_java') ~= 2; ...
      addpath(['midi_lib']); end

% The main function takes a MIDI file and an audio file as input
% The third arg is a flag to plot figures, and the fourth is a stem
% for the output midis.
midifile = 'SoS.mid';
wavfile = '06_Sultans_of_Swing.mp3';
doplot = 1;
stem = 'tmp';
alignmidi(midifile,wavfile,doplot,stem);
% With doplot set, the routine already plotted the results. zoom a little
subplot(122)
axis([0 300 0 300])

% A midi file warped to match audio is written to <stem>-mix.mid
outmidi = [stem, '-mix.mid'];
% Let's read it in, and play in stereo
sr = 11025;
domono = 1;
[da,sr] = audioread(wavfile, sr, domono);
[dm,sr] = midireadasaudio(outmidi, sr, domono);
soundsc([da(1:20*sr), dm(1:20*sr)], sr);

%% Comparing Outputs
%
% We can plot log-frequency spectrograms to show the effect and
% quality of the alignment.  Note that the aligned MIDI is not just
% shifted in time, but also scaled appropriately to keep every beat
% time lined up.

% Read original audio
[da,sr] = audioread(wavfile,11025,1);
% Read unaligned MIDI as audio
[dm,sr] = midireadasaudio(midifile,11025,1);
% .. and aligned version
[dma,sr] = midireadasaudio(outmidi,11025,1);
% Now plot
subplot(311)
logfsgram(da,512,sr)
caxis([-40 40])
axis([0 10 0.5 73.5])
title('Original audio')
subplot(312)
logfsgram(dma,512,sr)
caxis([-40 40])
axis([0 10 0.5 73.5])
title('Aligned MIDI')
subplot(313)
logfsgram(dm,512,sr)
caxis([-40 40])
axis([2.45+[0 10] 0.5 73.5])
title('Original MIDI (offset 2.45 sec to line up start)')

%% Parameters
%
% <alignmidi.m alignmidi> accepts a number of different parameters as fields of a 
% "params" structure, passed as the fifth argument.  They are:

%        params.tightness = 800 - beat tracking rigidity
%        params.voxchan = 4     - MIDI channel to isolate as vox
%        params.usedtw = 0      - 1 = use DTW; 0 = use Viterbi
%        params.gulley = 0.33   - proportion of edges where DTW can start/end
%        params.horizwt = 0.5   - penalty factor for horiz/vert DTW steps
%        params.priorwidth = 0.33 - spread of prior on Viterbi initial state
%        params.txwidth = 2.0   - viterbi transition with scale
%        params.txfloor = 0.001 - worst-case viterbi jump probability

%% Installation
% 
% This package has been compiled for several targets 
% using the Matlab compiler.  You will also need 
% to download and install the Matlab Compiler Runtime (MCR) Installer. 
% Please see the table below:
%
% <html>
% <table border=1>
% <tr><th>Architecture</th><th>Compiled package</th><th>MCR Installer</th></tr>
% <tr><td>64 bit Linux</td>
% <td><a href="alignmidi_GLNXA64.zip">alignmidi_GLNXA64.zip</a></td>
% <td><a href="http://www.ee.columbia.edu/~dpwe/tmp/MCRInstaller_glnxa64.bin">Linux 64 bit MCR Installer</a></td></tr>
% <tr><td>64 bit MacOS</td>
% <td><a href="alignmidi_MACI64.zip">alignmidi_MACI64.zip</a></td>
% <td><a href="http://www.ee.columbia.edu/~dpwe/tmp/MCRInstaller.dmg">MACI64 MCR Installer</a></td></tr>
% </table></html>
% 
% The original Matlab code used to build this compiled target is 
% available at <http://www.ee.columbia.edu/~dpwe/resources/matlab/alignmidi>
%
% All sources are in the package <alignmidi-v@VER@.zip>.
%
% Feel free to contact me with any problems.

%% Notes
%
% The included function <audioread.m audioread> is able to read a
% wide range of sound file types, but relies on a number of other
% packages and/or support functions being installed.  Most obscure
% of these is  ReadSound, a MEX wrapper I wrote for the dpwelib
% sound file interface.  See the 
% <http://labrosa.ee.columbia.edu/matlab/audioread/ audioread homepage>
% for more details.
%
% To link to the Java code for reading and writing MIDI files, you
% need to "edit classpath.txt" within Matlab, as described in 
% <http://www.ee.columbia.edu/~csmit/karaoke_midi.html Christine's
% instructions>.  And then you'll have to restart Matlab.  
% Else you'll get an error like "Undefined variable
% "PianoRollViewParser" or class "PianoRollViewParser.parse"."
%
% This code needs to convert MIDI to audio; this is accomplished by 
% <midi2wav.m>.  What this function does depends on the
% architecture: on my Mac, it runs the AppleScript script
% <midi2aiff.scpt>, which uses QuickTime Player 7.app to import
% MIDI and export a waveform.  If you don't have QTP7 installed, it
% won't work.  On Linux, it uses the open-source MIDI synthesizer 
% <http://timidity.sourceforge.net/ timidity>, which works well 
% *as long as you change the default sound font*
% to be fluidr3_gm.cfg 
% instead of freepats.cfg, by editing /etc/timidity/timidity.cfg.
% The bass drum in the freepats sound font is inexplicably awful.
%
% This code supercedes an earlier effort, 
% <http://labrosa.ee.columbia.edu/~dpwe/resources/matlab/alignmidiwav/ alignmidiwav> .

%% Referencing
% If you use this work in a publication, I would be grateful 
% if you referenced this page as follows:
%
% D. P. W. Ellis (2013).  "Aligning MIDI files to music audio", web resource.
% http://www.ee.columbia.edu/~dpwe/resources/matlab/alignmidi/

%% Acknowledgment
% This project was supported in part by the NSF under 
% grant IIS-1117015. Any opinions, findings and conclusions 
% or recommendations expressed in this material are those of the 
% authors and do not necessarily reflect the views of the Sponsors.

%% Changelog

% v0.04  2014-03-14 - Default alignment is now viterbi (which
%                     permits arbitrarily large jumps both
%                     backwards and forwards if the likelihood
%                     gains are large enough) instead of the
%                     strictly monotonic paths of DTW.  But DTW is
%                     much faster, and can still be selected with 
%                     params.usedtw=1
%
% v0.03  2014-03-12 - Major fixes to transposition estimation, DP penalty
%                     scoring, correct weighting of paths in
%                     gulleys, contrast binarizing of CQ spectra, 
%                     and promoting consistency between
%                     tempos in MIDI and audio. 
%
% v0.02  2013-07-19 - beat2 modified to use full bandwidth when
%                     estimating onset strength (was losing hi-hats
%                     by downsampling to 8kHz).  Small changes to
%                     beat timing in beat2.
%
% v0.01  2013-07-17 - Initial release

% Last updated: $Date: 2013/08/24 13:04:10 $
% Dan Ellis <dpwe@ee.columbia.edu>
