function data = wfft(data)
% WFFT One-sided Fourier transform of time domain data
%
% WFFT computes the one-sided Fourier transform of time domain data similar to
% the output produced by WCONDITION and WSCANCONDITION.  The resulting one-sided
% frequency domain data are assumed to extend from zero frequency to the Nyquist
% frequency.  As a result, time domain data of length N corresponds to frequency
% domain data of length N / 2 + 1.
%
% WFFT can be used as a placebo replacement for WCONDITION.
%
% usage:
%
% frequencyDomainData = wfft(timeDomainData);
%
%   timeDomainData        cell array of time domain data
%
%   frequencyDomainData   cell array of one-sided frequency domain data
%
% See also WIFFT, WCONDITION, and WSCANCONDITION.

% Shourov K. Chatterji
% shourov@ligo.mit.edu

% $Id: wfft.m 986 2008-08-14 21:04:56Z lstein $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for sufficient command line arguments
error(nargchk(1, 1, nargin));

% force cell array of data
data = wmat2cell(data);

% force one dimensional cell array
data = data(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(data);

% force row vectors
for channelNumber = 1 : numberOfChannels,
  data{channelNumber} = data{channelNumber}(:).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            fast fourier transform                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inverse fourier transform
for channelNumber = 1 : numberOfChannels,
  data{channelNumber} = fft(data{channelNumber});
  data{channelNumber} = ...
      data{channelNumber}(1 : length(data{channelNumber}) / 2 + 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
