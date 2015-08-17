function [Outfilt] = mhighpass(gain,samplFreq,order,freqCutOff)
%
% MHIGHPASS return a high pass filter suitable for use with "mfilter"
% mex file. A Butterworth filter is used.
% 
% Usage: Outfilt = mhighpass(gain,samplFreq,order,freqCutOff)
% 
% gain          : High-frequency gain.
% samplFreq     : Sampling frequency.
% order         : Order of the filter.
% freqCutOff    : Cut off frequency in Hz (0 < freqCutOff < samplFreq/2).
%     
% Outfilt       : A filter structure for use in mfilter.
%     
% M. Hewitson
%
% $Id: mhighpass.m 162 2009-08-14 20:11:29Z ajith $

if(freqCutOff > samplFreq/2)
  error('### freqCutOff must be < samplFreq/2');
end

Outfilt.name           = 'outfilt';
Outfilt.fs             = samplFreq;
Outfilt.ntaps          = order+1;
[Outfilt.a, Outfilt.b] = butter(order, 2*freqCutOff/samplFreq, 'high');

Outfilt.gain           = gain;
Outfilt.histin         = zeros(1,Outfilt.ntaps-1);       % initialise input history
Outfilt.histout        = zeros(1,Outfilt.ntaps-1);       % initialise output history



