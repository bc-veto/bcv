function varargout = wcondition(rawData, tiling, doubleWhiten)
% WCONDITION High pass filter and whiten time series data
%
% WCONDITION high pass filters and whitens time series data prior to analysis
% by the Q transform.  The data is first zero-phase high pass filtered at the
% minimum analysis frequency specified in the tiling structure.  The resulting
% high pass filtered data is then whitened by zero-phase linear prediction at
% a frequency resolution equal to the minimum analysis bandwidth requested in
% the tiling structure.  Note that the resulting whitened data is returned as
% a frequency series, not a time series.  In addition, the returned frequency
% series extends from zero frequency to the Nyquist frequency.  As a result,
% they are of length N / 2 + 1, where N is the length of the individual input
% time series.  The returned frequency-domain whitened data is normalized to
% unity amplitude spectral density.
%
% To enable recovery of calibrated amplitude and phase information, WCONDITION
% also returns the effective frequency domain coefficients of the combined high
% pass and whitening filters for each channel.
%
% WCONDITION also permits double whitening of the data to support true matched
% filtering.  Regardless of whether double whitening is requested or not, the
% returned frequency-domain filter coefficients always correspond to the single
% whitened case.
%
% usage:
%
%   [conditionedData, coefficients] = wcondition(rawData, tiling, doubleWhiten);
%
%   rawData               cell array of input time series
%   tiling                discrete Q transform tiling structure from WTILE
%
%   conditionedData       cell array of conditioned output frequency series
%   coefficients          cell array of frequency domain filter coefficients
%   doubleWhiten          boolean flag to enable double whitening
%
% There is also an alternative output syntax, which provides access to the
% intermediate raw and high pass filtered data for use by WSCAN.
%
%   [rawData, highPassedData, whitenedData, coefficients] = ...
%       wcondition(rawData, tiling, doubleWhiten);
%
%   rawData               cell array of unconditioned frequency series
%   highPassedData        cell array of high pass filtered frequency series
%   whitenedData          cell array of whitened frequency series data
%   coefficients          cell array of frequency domain filter coefficients
%
% See also WTILE, WTRANSFORM, WTHRESHOLD, WSELECT, WSEARCH, WEVENT, WSCAN,
% and WEXAMPLE.

% Shourov K. Chatterji <shourov@ligo.mit.edu>
% Antony Searle <antony.searle@anu.edu.au>
% Jameson Rollins <jrollins@phys.columbia.edu>

% $Id: wcondition.m 1857 2009-05-25 14:50:43Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 3, nargin));

% apply default arguments
if (nargin < 3) || isempty(doubleWhiten),
  doubleWhiten = 0;
end

% verify correct number of output arguments
error(nargchk(1, 4, nargout));

% determine necessary output arguments
switch nargout,
  case 1,
    returnIntermediateData = 0;
    returnCoefficients = 0;
  case 2,
    returnIntermediateData = 0;
    returnCoefficients = 1;
  case 3,
    returnIntermediateData = 1;
    returnCoefficients = 0;
  case 4,
    returnIntermediateData = 1;
    returnCoefficients = 1;
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

% validate tiling structure
if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
  error('input argument is not a discrete Q transform tiling structure');
end

% determine required data lengths
dataLength = tiling.sampleFrequency * tiling.duration;
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
lpefOrder = ceil(tiling.sampleFrequency * tiling.whiteningDuration);

% linear predictor error filter length
filterLength = lpefOrder + 1;

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

% initialize cell array of whitened data vectors
whitenedData = cell(numberOfChannels, 1);

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          fast fourier transform                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % fourier transform high pass filtered data
  highPassedData{channelNumber} = fft(highPassedData{channelNumber}) / sqrt(length(highPassedData{channelNumber}));

  % fourier transform raw data
  if returnIntermediateData,
    rawData{channelNumber} = fft(rawData{channelNumber}) / sqrt(length(rawData{channelNumber}));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      linear predictor error filter                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % if linear predictive whitening is requested,
  if tiling.whiteningDuration > 0,

    % estimate frequency-domain power spectrum
    powerSpectrum = real(highPassedData{channelNumber}).^2 + ...
                    imag(highPassedData{channelNumber}).^2;

    % estimate time-domain auto-correlation
    autoCorrelation = real(ifft(powerSpectrum));
    
    % solve for linear predictor error filter coefficients
    lpefCoefficients = levinson(autoCorrelation(1 : filterLength), lpefOrder);
    
    % create zero-phase frequency-domain whitening filter
    lpefCoefficients = [lpefCoefficients zeros(1, dataLength - filterLength)];
    lpefCoefficients = abs(fft(lpefCoefficients));

    % extract one-sided frequency-domain whitening filter coefficients
    lpefCoefficients = lpefCoefficients(1 : halfDataLength);

    % extract one-sided frequency-domain high pass filtered data
    highPassedData{channelNumber} = ...
        highPassedData{channelNumber}(1 : halfDataLength);

    % extract one-sided frequency-domain raw data
    if returnIntermediateData,
      rawData{channelNumber} = rawData{channelNumber}(1 : halfDataLength);
    end

    if tiling.highPassCutoff > 0
        lpefCoefficients(1:highPassCutoffIndex) = 0;
    end
        
    % include whitening filter in conditioning coefficients
    if returnCoefficients,
      coefficients{channelNumber} = ...
          coefficients{channelNumber} .* lpefCoefficients;
    end
  
    % apply whitening filter
    whitenedData{channelNumber} = ...
        lpefCoefficients .* highPassedData{channelNumber};
    
    % mean amplitude spectral density
    % normalization = sqrt(dataLength * nyquistFrequency * halfDataLength / ...
    %                      sum(abs(whitenedData{channelNumber}).^2));
    if tiling.highPassCutoff > 0
        normalization = 1 ./ sqrt(mean(abs(whitenedData{channelNumber}((highPassCutoffIndex + 1):end)).^2));
    else
        normalization = 1 ./ sqrt(mean(abs(whitenedData{channelNumber}).^2));
    end
    % under this normalization convention, <|x|^2> = 1, <re(x)^2> = 0.707 

    % renormalize whitened data to unity amplitude spectral density
    whitenedData{channelNumber} = ...
        normalization * whitenedData{channelNumber};

    % include normalization factor in whitening filter coefficients
    if returnCoefficients,
      coefficients{channelNumber} = ...
          coefficients{channelNumber} .* normalization;
    end

    % reapply whitening filter and normalization if double whitening requested
    if doubleWhiten,
      whitenedData{channelNumber} = ...
          normalization * lpefCoefficients .* whitenedData{channelNumber};
    end

  % if linear predictive whitening is not requested,
  else

    % extract one-sided frequency-domain high pass filtered data
    highPassedData{channelNumber} = ...
        highPassedData{channelNumber}(1 : halfDataLength);

    % extract one-sided frequency-domain raw data
    if returnIntermediateData,
      rawData{channelNumber} = ...
          rawData{channelNumber}(1 : halfDataLength);
    end

    % do nothing
    whitenedData{channelNumber} = highPassedData{channelNumber};
  
  % end test for linear predictive whitening
  end

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
    varargout{1} = whitenedData;
  case 2,
    varargout{1} = whitenedData;
    varargout{2} = coefficients;
  case 3,
    varargout{1} = rawData;
    varargout{2} = highPassedData;
    varargout{3} = whitenedData;
  case 4,
    varargout{1} = rawData;
    varargout{2} = highPassedData;
    varargout{3} = whitenedData;
    varargout{4} = coefficients;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
