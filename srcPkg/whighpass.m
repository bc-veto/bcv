function varargout = whighpass(rawData, tiling)
% WHIGHPASS High pass filter time series data
%
% WHIGHPASS high pass filters time series data. The data is zero-phase high pass 
% filtered at the minimum analysis frequency specified in the tiling structure.  
%
% usage:
%
%   [highPassedData, coefficients] = whighpass(rawData, tiling);
%
%   rawData               cell array of input time series
%   tiling                discrete Q transform tiling structure from WTILE
%
%   coefficients          cell array of frequency domain filter coefficients
%
%
% Shourov K. Chatterji <shourov@ligo.mit.edu>
% Antony Searle <antony.searle@anu.edu.au>
% Jameson Rollins <jrollins@phys.columbia.edu>
%
% $Id: whighpass.m 162 2009-08-14 20:11:29Z ajith $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 2, nargin));

% verify correct number of output arguments
error(nargchk(1, 2, nargout));

% determine necessary output arguments
switch nargout,
  case 1,
    returnCoefficients = 0;
  case 2,
    returnCoefficients = 1;
  otherwise,
      error('incorrect number of arguments');
end

% force cell arrays
if ~iscell(rawData),
  rawData = mat2cell(rawData, size(rawData, 1), size(rawData, 2));
end

% force one dimensional cell arrays
rawData = rawData(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(rawData);

% determine required data lengths
%dataLength = tiling.sampleFrequency * tiling.duration;
dataLength = length(rawData{1});
tiling.duration = dataLength/tiling.sampleFrequency;

halfDataLength = dataLength / 2 + 1;

% validate data length and force row vectors
for channelNumber = 1 : numberOfChannels,
  rawData{channelNumber} = rawData{channelNumber}(:).';
  if length(rawData{channelNumber}) ~= dataLength,
    error('data length not consistent with tiling');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                design filters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nyquist frequency
nyquistFrequency = tiling.sampleFrequency / 2;

% linear predictor error filter order
%lpefOrder = ceil(tiling.sampleFrequency * tiling.whiteningDuration);
lpefOrder = 0;

% if high pass filtering is requested,
if tiling.highPassCutoff > 0,
    
  % high pass filter order
  hpfOrder = 12;

  % design high pass filter
  [hpfZeros, hpfPoles, hpfGain] = ...
      butter(hpfOrder, tiling.highPassCutoff / nyquistFrequency, 'high');
  hpfSOS = zp2sos(hpfZeros, hpfPoles, hpfGain);

  % magnitude response of high pass filter
  minimumFrequencyStep = 1 / tiling.duration;
  frequencies = 0 : minimumFrequencyStep : nyquistFrequency;
  hpfArgument = (frequencies / tiling.highPassCutoff).^(2 * hpfOrder);
  hpfResponse = hpfArgument ./ (1 + hpfArgument);

  highPassCutoffIndex = ceil(tiling.highPassCutoff / minimumFrequencyStep);
    
% end test for high pass filtering
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            initialize cell arrays                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cell array of high pass filtered data vectors
highPassedData = cell(numberOfChannels, 1);

% initialize cell array of conditioning filter coefficients
if returnCoefficients,
  coefficients = cell(numberOfChannels, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                   initialize conditioning coefficients                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % initialize conditioning filter coefficients
  if returnCoefficients,
    coefficients{channelNumber} = ones(1, halfDataLength);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                             high pass filter                               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % if high pass filtering is requested,
  if tiling.highPassCutoff > 0,

    % apply high pass filter
    highPassedData{channelNumber} = sosfiltfilt(hpfSOS, rawData{channelNumber});

    % include high pass filter in conditioning coefficients
    if returnCoefficients,
      coefficients{channelNumber} = coefficients{channelNumber} .* hpfResponse;
    end
    
  % if high pass filtering is not requested,
  else,

    % do nothing
    highPassedData{channelNumber} = rawData{channelNumber};

  % end test for high pass filtering
  end

  % supress high pass filter transients
  highPassedData{channelNumber}(1 : lpefOrder) = ...
      zeros(1, lpefOrder);
  highPassedData{channelNumber}(dataLength - lpefOrder + 1 : dataLength) = ...
      zeros(1, lpefOrder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          construct output arguments                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct output arguments
switch nargout,
  case 1,
    varargout{1} = highPassedData;
  case 2,
    varargout{1} = highPassedData;
    varargout{2} = coefficients;
  otherwise, 
    error('Incorrect number of outputs');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;