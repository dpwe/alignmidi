function [Y,FS,NBITS,OPTS] = mp3read(FILE,N,MONO,DOWNSAMP,DELAY)
% MP3READ   Read MP3 audio file via use of external binaries.
%   Y = MP3READ(FILE) reads an mp3-encoded audio file into the
%     vector Y just like wavread reads a wav-encoded file (one channel 
%     per column).  Extension ".mp3" is added if FILE has none.
%     Also accepts other formats of wavread, such as
%   Y = MP3READ(FILE,N) to read just the first N sample frames (N
%     scalar), or the frames from N(1) to N(2) if N is a two-element vector.  
%   Y = MP3READ(FILE,FMT) or Y = mp3read(FILE,N,FMT) 
%     with FMT as 'native' returns int16 samples instead of doubles; 
%     FMT can be 'double' for default behavior (to exactly mirror the
%     syntax of wavread).
%
%   [Y,FS,NBITS,OPTS] = MP3READ(FILE...) returns extra information:
%     FS is the sampling rate,  NBITS is the bit depth (always 16), 
%     OPTS.fmt is a format info string; OPTS has multiple other
%     fields, see WAVREAD.
%
%   SIZ = MP3READ(FILE,'size') returns the size of the audio data contained
%     in the file in place of the actual audio data, returning the
%     2-element vector SIZ=[samples channels].
%
%   [Y...] = MP3READ(FILE,N,MONO,DOWNSAMP,DELAY) extends the
%     WAVREAD syntax to allow access to special features of the
%     mpg123 engine:  MONO = 1 forces output to be mono (by
%     averaging stereo channels); DOWNSAMP = 2 or 4 downsamples by 
%     a factor of 2 or 4 (thus FS returns as 22050 or 11025
%     respectively for a 44 kHz mp3 file); 
%     To accommodate a bug in mpg123-0.59, DELAY controls how many
%     "warm up" samples to drop at the start of the file; the
%     default value of 2257 makes an mp3write/mp3read loop for a 44
%     kHz mp3 file be as close as possible to being temporally
%     aligned; specify as 0 to prevent discard of initial samples.
%     For later versions of mpg123 (e.g. 1.9.0) this is not needed; 
%     a flag in mp3read.m makes the default DELAY zero in this case.
%
%   [Y...] = MP3READ(URL...)  uses the built-in network
%     functionality of mpg123 to read an MP3 file across the
%     network.  URL must be of the form 'http://...' or
%     'ftp://...'.  'size' and OPTS are not available in this mode.
%
%   Example:
%   To read an mp3 file as doubles at its original width and sampling rate:
%     [Y,FS] = mp3read('piano.mp3');
%   To read the first 1 second of the same file, downsampled by a
%   factor of 4, cast to mono, using the default filename
%   extension:
%     [Y,FS4] = mp3read('piano', FS/4, 1, 4);
%
%   Note: Because the mp3 format encodes samples in blocks of 26 ms (at
%   44 kHz), and because of the "warm up" period of the encoder,
%   the file length may not be exactly what you expect, depending 
%   on your version of mpg123 (recent versions fix warmup).
%
%   Note: requires external binaries mpg123 and mp3info; you
%   can find binaries for several platforms at:
%     http://labrosa.ee.columbia.edu/matlab/mp3read.html
%
%   See also mp3write, wavread.

% $Header: /Users/dpwe/matlab/columbiafns/RCS/mp3read.m,v 1.7 2010/04/09 18:13:00 dpwe Exp dpwe $

% 2003-07-20 dpwe@ee.columbia.edu  This version calls mpg123.
% 2004-08-31 Fixed to read whole files correctly
% 2004-09-08 Uses mp3info to get info about mp3 files too
% 2004-09-18 Reports all mp3info fields in OPTS.fmt; handles MPG2LSF sizes
%            + added MONO, DOWNSAMP flags, changed default behavior.
% 2005-09-28 Fixed bug reading full-rate stereo as 1ch (thx bjoerns@vjk.dk)
% 2006-09-17 Chop off initial 2257 sample delay (for 44.1 kHz mp3)
%            so read-write loop doesn't get progressively delayed.
%            You can suppress this with a 5th argument of 0.
% 2007-02-04 Added support for FMT argument to match wavread
%            Added automatic selection of binary etc. to allow it
%            to work cross-platform without editing prior to
%            submitting to Matlab File Exchange
% 2007-07-23 Tweaks to 'size' mode so it exactly agrees with read data.
% 2009-03-15 Added fixes so 'http://...' file URLs will work.
% 2009-03-26 Added filename length check to http: test (thx fabricio guzman)

% find our baseline directory
persistent path mpg123 mp3info
if isempty(path)
  path = fileparts(which('mp3read'));
end

% %%%%% Directory for temporary file (if needed)
% % Try to read from environment, or use /tmp if it exists, or use CWD
tmpdir = getenv('TMPDIR');
if isempty(tmpdir) || exist(tmpdir,'file')==0
  tmpdir = '/tmp';
end
if exist(tmpdir,'file')==0
  tmpdir = '';
end
% ensure it exists
%if length(tmpdir) > 0 && exist(tmpdir,'file')==0
%  mkdir(tmpdir);
%end

%%%%%% Command to delete temporary file (if needed)
rmcmd = 'rm';

%%%%%% Location of the binaries - attempt to choose automatically
%%%%%% (or edit to be hard-coded for your installation)
ext = lower(computer);
if ispc
  ext = 'exe';
  rmcmd = 'del';
end
% mpg123-0.59 inserts silence at the start of decoded files, which
% we compensate.  However, this is fixed in mpg123-1.9.0, so 
% make this flag 1 only if you have mpg123-0.5.9
MPG123059 = 0;
%mpg123 = fullfile(path,['mpg123.',ext]);
%mp3info = fullfile(path,['mp3info.',ext]);
if isempty(mpg123)
  [r,mpg123] = system('which mpg123');
  if r ~= 0; error(mpg123); end
  [r,mp3info] = system('which mp3info');
  if r ~= 0; error(mp3info); end
  % strip any returns
  mpg123 = mpg123(double(mpg123)>31);
  mp3info = mp3info(double(mp3info)>31);
end
  
%%%%% Check for network mode
if length(FILE) > 6 && (strcmp(lower(FILE(1:7)),'http://') == 1 ...
      || strcmp(lower(FILE(1:6)),'ftp://'))
  % mp3info not available over network
  OVERNET = 1;
else
  OVERNET = 0;
end


%%%%% Process input arguments
if nargin < 2
  N = 0;
end

% Check for FMT spec (per wavread)
FMT = 'double';
if ischar(N)
  FMT = lower(N);
  N = 0;
end

if length(N) == 1
  % Specified N was upper limit
  N = [1 N];
end
if nargin < 3
  forcemono = 0;
else
  % Check for 3rd arg as FMT
  if ischar(MONO)
    FMT = lower(MONO);
    MONO = 0;
  end
  forcemono = (MONO ~= 0);
end
if nargin < 4
  downsamp = 1;
else
  downsamp = DOWNSAMP;
end
if downsamp ~= 1 && downsamp ~= 2 && downsamp ~= 4
  error('DOWNSAMP can only be 1, 2, or 4');
end

% process DELAY option (nargin 5) after we've read the SR

if strcmp(FMT,'native') == 0 && strcmp(FMT,'double') == 0 && ...
      strcmp(FMT,'size') == 0
  error(['FMT must be ''native'' or ''double'' (or ''size''), not ''',FMT,'''']);
end


%%%%%% Constants
NBITS=16;

%%%%% add extension if none (like wavread)
[path,file,ext] = fileparts(FILE);
if isempty(ext)
  FILE = [FILE, '.mp3'];
end

%%%%% maybe expand ~ %%%%%%
if FILE(1) == '~'
  FILE = [getenv('HOME'),FILE(2:end)];
end

%%%%% escape weird shell characters %%%%%
badchrs = '"\!`';
ix = find(max(repmat(FILE,length(badchrs),1)==repmat(badchrs',1,length(FILE))));
offset = 0;
for i = ix
  % Move string from "bad" character to end down one space
  FILE = [FILE,' '];
  FILE((i+offset)+1:end) = FILE((i+offset):end-1);
  % insert a backslash before the "bad" character
  FILE(i+offset) = '\';
  % remember that all the subsequent characters are now one 
  % character further down
  offset = offset+1;
end


if ~OVERNET
  %%%%%% Probe file to find format, size, etc. using "mp3info" utility
  cmd = ['"',mp3info, '" -r m -p "%Q %u %b %r %v * %C %e %E %L %O %o %p" "', FILE,'"'];
  % Q = samprate, u = #frames, b = #badframes (needed to get right answer from %u) 
  % r = bitrate, v = mpeg version (1/2/2.5)
  % C = Copyright, e = emph, E = CRC, L = layer, O = orig, o = mono, p = pad
  w = mysystem(cmd);
  % Break into numerical and ascii parts by finding the delimiter we put in
  starpos = findstr(w,'*');
  nums = str2num(w(1:(starpos - 2)));
  strs = tokenize(w((starpos+2):end));

  SR = nums(1);
  nframes = nums(2);
  nchans = 2 - strcmp(strs{6}, 'mono');
  layer = length(strs{4});
  bitrate = nums(4)*1000;
  mpgv = nums(5);
  % Figure samples per frame, after
  % http://board.mp3-tech.org/view.php3?bn=agora_mp3techorg&key=1019510889
  if layer == 1
    smpspfrm = 384;
  elseif SR < 32000 && layer ==3
    smpspfrm = 576;
    if mpgv == 1
      error('SR < 32000 but mpeg version = 1');
    end
  else
    smpspfrm = 1152;
  end

  OPTS.fmt.mpgBitrate = bitrate;
  OPTS.fmt.mpgVersion = mpgv;
  % fields from wavread's OPTS
  OPTS.fmt.nAvgBytesPerSec = bitrate/8;
  OPTS.fmt.nSamplesPerSec = SR;
  OPTS.fmt.nChannels = nchans;
  OPTS.fmt.nBlockAlign = smpspfrm/SR*bitrate/8;
  OPTS.fmt.nBitsPerSample = NBITS;
  OPTS.fmt.mpgNFrames = nframes;
  OPTS.fmt.mpgCopyright = strs{1};
  OPTS.fmt.mpgEmphasis = strs{2};
  OPTS.fmt.mpgCRC = strs{3};
  OPTS.fmt.mpgLayer = strs{4};
  OPTS.fmt.mpgOriginal = strs{5};
  OPTS.fmt.mpgChanmode = strs{6};
  OPTS.fmt.mpgPad = strs{7};
  OPTS.fmt.mpgSampsPerFrame = smpspfrm;
else
  % OVERNET mode
  OPTS = [];
  % guesses
  smpspfrm = 1152;
  SR = 44100;
  nchans = 2;
  nframes = 0;
end
  
if SR == 16000 && downsamp == 4
  error('mpg123 will not downsample 16 kHz files by 4 (only 2)');
end

% from libmpg123/frame.h
GAPLESS_DELAY = 529;

% process or set delay
if nargin < 5

  if MPG123059
    mpg123delay44kHz = 2257;  % empirical delay of lame/mpg123 loop
    mpg123delay16kHz = 1105;  % empirical delay of lame/mpg123 loop
                              % for 16 kHz sampling - one 1152
                              % sample frame less??
    if SR == 16000
      rawdelay = mpg123delay16kHz;
    else
      rawdelay = mpg123delay44kHz;  % until we know better
    end
    delay = round(rawdelay/downsamp);
  else
    % seems like predelay is fixed in mpg123-1.9.0
    delay = 0;
  end
else
  delay = DELAY;
end

if downsamp == 1
  downsampstr = '';
else
  downsampstr = [' -',num2str(downsamp)];
end
FS = SR/downsamp;

if forcemono == 1
  nchans = 1;
  chansstr = ' -m';
else
  chansstr = '';
end

% Size-reading version
if strcmp(FMT,'size') == 1
  if OVERNET
    % no info over net
    Y = [-1 0];
  else
    if MPG123059
      Y = [floor(smpspfrm*nframes/downsamp)-delay, nchans];
    else
      Y = [floor(smpspfrm*nframes/downsamp)-GAPLESS_DELAY, nchans];
    end
  end
else

  % Temporary file to use
  tmpfile = fullfile(tmpdir, ['tmp',num2str(round(1000*rand(1))),'.wav']);

  skipx = 0;
  skipblks = 0;
  skipstr = '';
  sttfrm = N(1)-1;

  % chop off transcoding delay?
  %sttfrm = sttfrm + delay;  % empirically measured
  % no, we want to *decode* those samples, then drop them
  % so delay gets added to skipx instead
  
  if sttfrm > 0
    skipblks = floor(sttfrm*downsamp/smpspfrm);
    skipx = sttfrm - (skipblks*smpspfrm/downsamp);
    skipstr = [' -k ', num2str(skipblks)];
  end
  skipx = skipx + delay;
  
  lenstr = '';
  endfrm = -1;
  decblk = 0;
  if length(N) > 1 && N(2) > 0
    endfrm = N(2);
    if endfrm > sttfrm
      decblk = ceil((endfrm+delay)*downsamp/smpspfrm) - skipblks + 10;   
      % we read 10 extra blks (+10) to cover the case where up to 10 bad 
      % blocks are included in the part we are trying to read (it happened)
      lenstr = [' -n ', num2str(decblk)];
      % This generates a spurious "Warn: requested..." if reading right 
      % to the last sample by index (or bad blks), but no matter.
    end
 end

  % Run the decode
  cmd=['"',mpg123,'"', downsampstr, chansstr, skipstr, lenstr, ...
       ' -q -w "', tmpfile,'"  "',FILE,'"'];
  %w = 
%  disp(cmd);
  mysystem(cmd);

  % Load the data (may update FS if it was based on a guess previously)
  [Y,FS] = wavread(tmpfile);

%  % pad delay on to end, just in case
%  Y = [Y; zeros(delay,size(Y,2))];
%  % no, the saved file is just longer
  
  if decblk > 0 && length(Y) < decblk*smpspfrm/downsamp
    % This will happen if the selected block range includes >1 bad block
    disp(['Warn: requested ', num2str(decblk*smpspfrm/downsamp),' frames, returned ',num2str(length(Y))]);
  end
  
  % Delete tmp file
  mysystem([rmcmd,' "', tmpfile,'"']);
  
  % debug
%  disp(['sttfrm=',num2str(sttfrm),' endfrm=',num2str(endfrm),' skipx=',num2str(skipx),' delay=',num2str(delay),' len=',num2str(length(Y))]);
  
  % Select the desired part
  if skipx+endfrm-sttfrm > length(Y)
      endfrm = length(Y)+sttfrm-skipx;
  end

  if endfrm > sttfrm
    Y = Y(skipx+(1:(endfrm-sttfrm)),:);
  elseif skipx > 0
    Y = Y((skipx+1):end,:);
  end
  
  % Convert to int if format = 'native'
  if strcmp(FMT,'native')
    Y = int16((2^15)*Y);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0 
  error(['unable to execute ',cmd,' (',w(1:end-1),')']);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
% Debug
%disp([cmd,' -> ','*',w,'*']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = tokenize(s,t)
% Break space-separated string into cell array of strings.
% Optional second arg gives alternate separator (default ' ')
% 2004-09-18 dpwe@ee.columbia.edu
if nargin < 2;  t = ' '; end
a = [];
p = 1;
n = 1;
l = length(s);
nss = findstr([s(p:end),t],t);
for ns = nss
  % Skip initial spaces (separators)
  if ns == p
    p = p+1;
  else
    if p <= l
      a{n} = s(p:(ns-1));
      n = n+1;
      p = ns+1;
    end
  end
end
    
