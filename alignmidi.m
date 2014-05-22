function [p,q,SM,DAB,DMB,baa,bma] = alignmidi(MF,WF,do_plot,out,params)
% [p,q,S,D,M,ba,bm] = alignmidi(MF,WF,doplot,out,params)
%     Align a midi file to a wav file by resynthesizing, beat
%     tracking, and DP over logsgram from 41.2 to 2 kHz over
%     whitenened waveforms.
%     Inputs:
%      MF is the name of the MIDI file
%      WF is the name of the wav file.
%      doplot = 1 plots the DTW path
%      out is the stem for the outputs.
%      params holds various other parameters:
%        params.tightness = 800 - beat tracking rigidity
%        params.voxchan = 4     - MIDI channel to isolate as vox
%        params.usedtw = 0      - 1 = use DTW; 0 = use Viterbi
%        params.gulley = 0.33   - proportion of edges where DTW can start/end
%        params.horizwt = 0.5   - penalty factor for horiz/vert DTW steps
%        params.priorwidth = 0.33 - spread of prior on Viterbi initial state
%        params.txwidth = 2.0   - viterbi transition with scale
%        params.txfloor = 0.001 - worst-case viterbi jump probability
%     Outputs:
%      [p,q] are the path from DP
%      S is the similarity matrix.
%      D is the logf-spectrogram of the wav, M is for the midi.
%      ba, bm are the aligned, corresponding beat times in audio
%       and midi respectively.
%  2013-07-16 Dan Ellis dpwe@ee.columbia.edu

VERSION = 0.04;
DATE = 20140314;

if nargin < 3;  do_plot = 0; end
if nargin < 4;  out = MF; end
if nargin < 5;  params.foo = []; end

% minimal support for command-line usage
if ischar(do_plot);  do_plot = str2num(do_plot); end
% .. but not for params...

% tightness for beat2 - larger means more rigid tempo
if ~isfield(params, 'tightness'); params.tightness = 800; end
% MIDI channel to isolate for vocals - 4 by convention
if ~isfield(params, 'voxchan');   params.voxchan = 4; end
% Do we use DTW or Viterbi alignment?
if ~isfield(params, 'usedtw');    params.usedtw = 0; end
% Proportion of initial and final edges to accept for DTW
% (mismatched start/end)
if ~isfield(params, 'gulley');    params.gulley = 0.3; end
% Extra weighting for horizontal/vertical steps in DTW (to
% encourage diagonal)
if ~isfield(params, 'horizwt');   params.horizwt = 0.5; end
% Equivalent of gulley for Viterbi - how broad is initial state prior
if ~isfield(params, 'priorwidth');   params.priorwidth = 0.5; end
% Equivalent horizwt for Viterbi - how wide is laplacian on transition
if ~isfield(params, 'txwidth');   params.txwidth = 2.0; end
% Probability of jump to any remote spot for Viterbi
if ~isfield(params, 'txfloor');   params.txfloor = 0.1; end


% output file names
[p,n,e] = fileparts(out);
MFout = fullfile(p, [n, '-mix.mid']);
MFvox = fullfile(p, [n, '-vox.mid']);
MFins = fullfile(p, [n, '-ins.mid']);

% Read midi synthesized audio
sr = 11025;
dm = midireadasaudio(MF, sr, 1);

% Read actual audio
da = audioread(WF, sr, 1);

% Calculate beat times and beat-sync log-f sgram features
[DMB,bm,ff] = beatsynclogspec(dm, sr);
Mbpm = 60/median(diff(bm));
[DAB,ba,ff] = beatsynclogspec(da, sr);
Abpm = 60/median(diff(ba));
disp(['Tempo: MIDI=',num2str(Mbpm),' audio=',num2str(Abpm)]);
if abs(log(Mbpm/Abpm)) > (log(1.2))
  disp('**** WARNING: likely tempo mismatch - redoing MIDI');
  % Redo midi using audio tempo?
  [DMB,bm,ff] = beatsynclogspec(dm, sr, [Abpm, .5]);
  Mbpm = 60/median(diff(bm));
  disp(['Tempo: MIDI=',num2str(Mbpm),' audio=',num2str(Abpm)]);  
end


% Check for transposition
dropbass = 40;  % lowest channels are blurry & high energy - don't use
midisharpsemis = findtransposition(DAB, DMB, dropbass);
if midisharpsemis ~= 0
  disp(['Transposition detected: MIDI is ',num2str(midisharpsemis), ...
        ' semitones sharp']);
  DMB = rot(DMB', midisharpsemis)';
end

% further compress / normalize features
ampexp = 0.3;
DMB = DMB.^ampexp;
DAB = DAB.^ampexp;
%DMB = rownorm01(DMB);
%DAB = rownorm01(DAB);

% similarity matrix on normalized features
nfcw = 51;
%SM =  1 - simmx(rownorm01(normftrcols(DAB,nfcw)), ...
%                rownorm01(normftrcols(DMB,nfcw)));
K = [-1 -2 -1; 2 4 2; -1 -2 -1];
SM = simmx(double(conv2(rownorm01(normftrcols(DMB,nfcw)),K,'same')>0), ...
           double(conv2(rownorm01(normftrcols(DAB,nfcw)),K,'same')>0));


if params.usedtw
  disp('Using traditional DTW');
  % best path
  %[p,q] = dpfast(SM,[[1 1 1.0; 0 1 horizwt;1 0 horizwt]],0,gulley);
  %horizwt = 1.1; % added cost of horizontal moves
  %gulley = 0.2;  % proportion of final edges OK
  [p,q,C,phi,score] = dpmod(1-SM', params.horizwt, params.gulley);
  lp = length(p);
  l1090 = (round(0.1*lp)+1):round(0.9*lp);
  disp(['SD of 10..90% of path = ',num2str(std(p(l1090)-q(l1090)))]);
  %disp(['DP Best cost per pt = ',num2str(score/length(p))]);
  disp(['DP Best cost per pt 10..90 = ', ...
        num2str(mean(SM(sub2ind(size(SM),q(l1090),p(l1090)))))]);
else % USEVITERBI
  disp('Using viterbi alignment');
  nr = size(SM,1);
  pri = exp(-0.5*([1:nr]/(params.priorwidth*nr)).^2);
  rtx = max(params.txfloor, exp(-abs([-(nr-2):(nr)]'/params.txwidth)));
  [pp,tc,trb] = viterbi_implicittx(SM, pri, rtx);
  p = 1:size(SM,2);
  q = pp;
end

% aligned beat times
bma = bm(q);
baa = ba(p);

if do_plot
  % graphical debug
  subplot(321)
  cols = 1:200;
  imgsc(ba(cols),1:size(DAB,1),DAB(:,cols));
  yt = get(gca,'YTick');
  for i = 1:length(yt)
    ytl{i} = sprintf('%.0f',ff(yt(i)));
  end
  set(gca,'YTickLabel',ytl);
  title(['Wavfile ',WF,' beat-sync logfsgram'], 'interpreter','none');
  ylabel('freq / Hz');
  xlabel('time / sec');
  
  subplot(323)
  imgsc(bm(cols),1:size(DAB,1),DMB(:,cols));
  set(gca,'YTickLabel',ytl);
  title(['Midifile ',MF,' beat-sync logfsgram'], 'interpreter','none');
  ylabel('freq / Hz');
  xlabel('time / sec');

  subplot(122)
  imgsc(SM);
  hold on; plot(p,q,'r'); hold off
  title('Similarity matrix and best path');
  xlabel('time in audio / beats')
  ylabel('time in midi / beats');
  
  subplot(325)
  plot(p, q-p);
  xlabel('beats in ref');
  ylabel('delay of midi / beats');
  
  gcolor;
end

% Rewrite MIDI with fixed times
nmat = readmidi_java(MF,true);
% Apply transposition
nmat(:,4) = nmat(:,4) - midisharpsemis;
% durations become end times
nmat(:,7) = nmat(:,6) + nmat(:,7);
% map the times
nmat(:,[6 7]) = maptimes_interp(nmat(:,6:7),bma,baa);
% end times back to durations
nmat(:,7) = nmat(:,7) - nmat(:,6);
% read in and modify the pitch bend
pchanges = get_pitch_changes(MF);
pchanges(:,3) = maptimes_interp(pchanges(:,3),bma,baa);
% read in the program changes
changes = get_program_changes(MF);
% write out
if length(out)
  writemidi_seconds(nmat,MFout,changes(:,[5 4 1]),pchanges(:,[3 4 5 1]));
  %writemidi_seconds(nmat,MFout,changes(:,[5 4 1]));
  disp(['Saved to ', MFout]);

  % write out separated parts
  %vxch = 4;
  vxch = params.voxchan;
  vxnotes = find(nmat(:,3) == vxch);
  if length(vxnotes) == 0
    disp(['Vox channel ', num2str(vxch),' is empty']);
  else
    writemidi_seconds(nmat(vxnotes,:),MFvox,changes(:,[5 4 1]), ...
                      pchanges(find(pchanges(:,4)==vxch),[3 4 5 1]));
    disp(['Wrote ', MFvox]);

    writemidi_seconds(nmat(find(nmat(:,3)~=vxch),:),MFins,changes(:,[5 4 1]), ...
                    pchanges(find(pchanges(:,4)~=vxch),[3 4 5 1]));
    disp(['Wrote ', MFins]);
  end
end
