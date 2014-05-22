function AFout = midi2wav(MF,AF)
% AFout = midi2wav(MF,AF)
%   Convert a MIDI file into an audio file AF.
%   This version uses QuickTime Player 7 via Applescript on Mac
%   or Timidity on other platforms (Linux).
%   The WAV file is actually an AIFF file on Mac.
%   If AF is not specified, use a tmp file name, which is returned.
% 2013-07-16 Dan Ellis dpwe@ee.columbia.edu

comp = computer();
ismac = strcmp(comp(1:3), 'MAC');

if ismac
  wavext = '.aiff';
else
  wavext = '.wav';
end

if nargin < 2;   AF = [tempname(),wavext]; end

if length(AF) == 0
  % Derive names from midi file
  [p,n,e] = fileparts(MF);
  AF = fullfile(p, [n, wavext]);
end

if ismac % although... now we have timidity on mac (thanks, homebrew!)
  % Use QuickTime on Mac

  % Make sure output file is not there (else QT script bombs)
  if exist(AF, 'file')
    delete(AF);
  end

  % Run the applescript
  system(['osascript midi2aiff.scpt "',FNexpand(MF),'" "',FNexpand(AF),'"']);

else
  % Use Timidity on Linux

  system(['timidity -Ow -o ',AF,' ',MF]);
  % -Oa is AIFF output, -Ow is MSWAVE output        

end

if nargout > 0
  AFout = AF;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = FNexpand(F)
% convert a file name into a simple absolute path

if F(1) ~= '/'
  if F(1) == '~'
    G = [getenv('HOME'),F(2:end)];
  else
    G = fullfile(pwd(),F);
  end
else
  G = F;
end
