function trans = findtransposition(DA, DB, dropbass, do_plot)
% trans = findtransposition(DA, DB, dropbass, do_plot)
%   Check for transposition between two log-frequency spectrograms
%   DA and DB.
% 2013-07-17 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; dropbass = 0; end
if nargin < 4; do_plot = 0;  end

% Check for transposition
% skip the bass because it's blurry
%dropbass = 40;
% Cross-correlate the average log-freq spectra to find best pitch
% alignment 
xc = xcorr(rownorm01(mean(DB((dropbass+1):end,:),2)), ...
           rownorm01(mean(DA((dropbass+1):end,:),2)));
[vv,xcp] = max(xc);
xclen = size(DA,1)-dropbass-1;
trans = (xcp - xclen - 1);
if do_plot
  disp(['Transposition detected: MIDI is ',num2str(trans), ...
        ' semitones sharp']);
  subplot(311)
  plot([-xclen:xclen],xc)
  subplot(312)
  plot(1:xclen+1,rownorm01(mean(DMB((dropbass+1):end,:),2)), ...
       1:xclen+1,rownorm01(mean(DAB((dropbass+1):end,:),2)));
  subplot(325)
  imgsc(DA(:,101:200));
  subplot(326)
  imgsc(DB(:,101:200));
end
