function [data, injection] = qinject(noise, sampleFrequencies, startTime, ...
                                     tiling, waveform, signalToNoiseRatio, ...
                                     spectrumBand, spectrumResolution)
% QINJECT Inject simulated gravitational waves into detector noise
%
% QINJECT injects simulated gravitational waves with a variety of waveform
% into detector noise.  The specified waveform is injected at the specified
% signal to noise ratio, while other signal parameters are randomly selected
% to cover the signal space of the search.
%
% usage: [signal, injection] = qinject(noise, sampleFrequencies, startTime, ...
%                                      tiling, waveform, signalToNoiseRatio, ...
%                                      spectrumBand, spectrumResolution);
%
%   noise                cell array of detector noise data
%   sampleFrequencies    vector of sample frequency of detector noise data
%   startTime            GPS start time of detector noise data
%   tiling               discrete Q transform tiling structure from QTILE
%   waveform             string name of waveform to inject
%   signalToNoiseRatio   signal to noise ratio to inject
%   spectrumBand         frequency band for SNR estimation [Hz]
%   spectrumResolution   frequency resolution for SNR estimation [Hz]
%
%   data                 cell array of combined noise and injection data
%   injection            cell array containing Q transform injection structure
%
% The following waveform types are currently available
%
%   gaussian        simple gaussian pulse
%                     logarithmically distributed random duration
%                     uniformly distributed random sign
%   sinegaussian    gaussian enveloped sinusoid
%                     logarithmically distributed random quality factor
%                     logarithmically distributed random central frequency
%                     uniformly distributed random phase
%   noiseburst      bandlimited and time windowed white noise burst
%                     logarithmically distributed random central frequency
%                     logarithmically distributed random duration
%                     logarithmically distributed random bandwidth
%   inspiral        first order post newtonian inspiral
%                     uniformly distributed random mass
%                     uniformly distributed random phase
%   ringdown        exponentially damped sinusoid
%                     logarithmically distributed random quality factor
%                     logarithmically distributed random central frequency
%                     uniformly distributed random phase
%
% The specified frequency band for spectrum estimation should be a two
% component vector in the form [minimumFrequency maximumFrequency].  The
% amplitude of the injected signal to set to achieve the specified matched
% filter signal to noise ratio in this frequency band.  If spectrumBand or
% spectrumResolution are not specified, they default to [40 7000] and 1 Hz
% respectively.
%
% QINJECT returns a cell array containing a single Q transform
% injection structures that describes the properties of the injected
% signal. The event structure contains the following fields.
%
%   time                 center time of injection [gps seconds]
%   frequency            center frequency of injection [Hz]
%   duration             duration of injection [seconds]
%   bandwidth            bandwidth of injection [seconds]
%   snr                  signal to noise ratio of injection []
%   amplitude            characteristic amplitude of injection [Hz^-1/2]
%   phase                phase of injection [radians]
%
% If no spectrumResolution is specified, a default value of 1 Hz is assumed.
%
% See also QTESTCLUSTER.

% NOTES:
%
% QINJECT currently assumes that the tiling parmaeters tiling.minimumFrequency
% and tiling.maximum are not 0 or Inf.

% Shourov K. Chatterji
% shourov@ligo.mit.edu
% 2006-Jul-16

% $Id: qinject.m 2263 2009-08-24 19:54:40Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(6, 8, nargin));

% apply default arguments
if nargin < 7,
  spectrumBand = [40 7000];
end
if nargin < 8,
  spectrumResolution = 1;
end

% if input noise is not a cell array,
if ~iscell(noise),

  % insert noise into a single cell
  noise = mat2cell(noise, size(noise, 1), size(noise, 2));

% otherwise, continue
end

% force one dimensional cell array
noise = noise(:);

% force row vector of sample frequencies
sampleFrequencies = sampleFrequencies(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(noise);

% validate number of sample frequencies
if length(sampleFrequencies) ~= numberOfChannels,
  error('number of sample frequencies and channels are inconsistent');
end

% validate tiling structure
if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
  error('input argument is not a discrete Q transform tiling structure');
end

% validate data durations
for channelNumber = 1 : numberOfChannels,
  if (tiling.duration ~= ...
      length(noise{channelNumber}) / sampleFrequencies(channelNumber)),
    error('detector noise data durations not consistent with tiling');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    timing                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum channel sample frequency
maximumSampleFrequency = max(sampleFrequencies);

% maximum detector noise data length
maximumDataLength = maximumSampleFrequency * tiling.duration;

% stop time of detector noise data
stopTime = startTime + tiling.duration;

% generate time vector at maximum channel sample frequency
time = startTime + (0 : 1 : maximumDataLength - 1) / maximumSampleFrequency;

% estimated duration of filter transients to avoid
transientDuration = 4 * tiling.minimumFrequency / tiling.maximumQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           switch on waveform type                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(waveform)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                 Gaussian                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'gaussian'

    % select random central time [seconds]
    maximumSignalDuration = 1;
    minimumCentralTime = startTime + transientDuration + maximumSignalDuration;
    maximumCentralTime = stopTime - transientDuration - maximumSignalDuration;
    centralTime = unifrnd(minimumCentralTime, maximumCentralTime);

    % select random duration [seconds]
    minimumDuration = 1 / (2 * sqrt(pi) * tiling.maximumFrequency);
    maximumDuration = 1 / (2 * sqrt(pi) * tiling.minimumFrequency);
    duration =  2.^(unifrnd(log2(minimumDuration), log2(maximumDuration)));
    
    % select random phase [radians]
    phase = pi * binornd(1, 0.5);

    % normalize for unity characteristic amplititude
    amplitude = (2 * pi * duration^2)^(-1/4);
    
    % gaussian signal
    gaussian = exp(-(time - centralTime).^2 / (4 * duration^2));

    % generate signal
    signal = real(exp(-sqrt(-1) * phase)) * amplitude * gaussian;

    % record injection properties
    properties.time = centralTime;
    properties.frequency = 0;
    properties.duration = 2 * sqrt(pi) * duration;
    properties.bandwidth = 1 / (2 * sqrt(pi) * duration);
    properties.amplitude = 1;
    properties.phase = phase;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                           sinusoidal Gaussian                              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'sinegaussian'

    % select random central time [seconds]
    maximumSignalDuration = 1;
    minimumCentralTime = startTime + transientDuration + maximumSignalDuration;
    maximumCentralTime = stopTime - transientDuration - maximumSignalDuration;
    centralTime = unifrnd(minimumCentralTime, maximumCentralTime);

    % select random central frequency [Hz]
    minimumCentralFrequency = tiling.minimumFrequency;
    maximumCentralFrequency = tiling.maximumFrequency;
    centralFrequency =  2.^(unifrnd(log2(minimumCentralFrequency), ...
                                    log2(maximumCentralFrequency)));

    % select random quality factor []
    minimumQualityFactor = tiling.minimumQ;
    maximumQualityFactor = tiling.maximumQ;
    qualityFactor =  2.^(unifrnd(log2(minimumQualityFactor), ...
                                 log2(maximumQualityFactor)));

    % select random phase [radians]
    minimumPhase = 0;
    maximumPhase = 2 * pi;
    phase = unifrnd(minimumPhase, maximumPhase);

    % normalize for unity characteristic amplitude
    amplitude = (32 * pi * centralFrequency^2 / qualityFactor^2)^(1/4) * ...
                (1 - cos(2 * phase) * exp(-qualityFactor^2 / 2))^(-1/2);

    % gaussian envelope
    envelope = exp(-4 * pi^2 * centralFrequency^2 * ...
                   (time - centralTime).^2 / qualityFactor^2);

    % sinusoidal oscillation
    sinusoid = sin(2 * pi * centralFrequency * (time - centralTime) + phase);

    % generate signal
    signal = amplitude * envelope .* sinusoid;

    % record injection properties
    properties.time = centralTime;
    properties.frequency = centralFrequency;
    properties.duration = qualityFactor / (2 * sqrt(pi) * centralFrequency);
    properties.bandwidth = 2 * sqrt(pi) * centralFrequency / qualityFactor;
    properties.amplitude = 1;
    properties.phase = phase;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                            white noise burst                               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'noiseburst'

    % select random central time
    maximumSignalDuration = 2;
    minimumCentralTime = startTime + transientDuration + maximumSignalDuration;
    maximumCentralTime = stopTime - transientDuration - maximumSignalDuration;
    centralTime = unifrnd(minimumCentralTime, maximumCentralTime);

    % select random central frequency
    minimumCentralFrequency = tiling.minimumFrequency;
    maximumCentralFrequency = tiling.maximumFrequency;
    centralFrequency =  2.^(unifrnd(log2(minimumCentralFrequency), ...
                                    log2(maximumCentralFrequency)));

    % constrain minimum time frequency area of burst
    minimumTimeFrequencyArea = 64 / (4 * pi);

    % constrain maximum duration of burst
    maximumDuration = 1 / (2 * sqrt(pi));

    % select random bandwidth
    minimumBandwidth = minimumTimeFrequencyArea / maximumDuration;
    maximumBandwidth = centralFrequency / sqrt(11);
    bandwidth =  2.^(unifrnd(log2(minimumBandwidth), log2(maximumBandwidth)));

    % select random duration
    minimumDuration = minimumTimeFrequencyArea / (4 * pi * bandwidth);
    duration =  2.^(unifrnd(log2(minimumDuration), log2(maximumDuration)));

    % handle odd length vector
    oddDataLength = 0;
    if mod(length(time), 2) ~= 0,
      time(end + 1) = time(end) + 1 / maximumSampleFrequency;
      oddDataLength = 1;
    end

    % generate frequency domain white noise
    signal = randn(1, length(time) / 2 + 1) + ...
             sqrt(-1) * randn(1, length(time) / 2 + 1);
    signal = [signal conj(fliplr(signal(2 : end - 1)))];

    % generate frequency vector
    frequency = 0 : 1 / tiling.duration : maximumSampleFrequency / 2;
    frequency = [frequency -fliplr(frequency(2 : end - 1))];

    % window in frequency
    envelope = exp(-(abs(frequency) - centralFrequency).^2 ./ (4 * bandwidth^2));
    signal = signal .* envelope;

    % convert to time domain
    signal = real(ifft(signal));

    % preserve vector orientation
    reshape(signal, size(time));

    % window in time
    envelope = exp(-(time - centralTime).^2 / (4 * duration^2));
    signal = envelope .* signal;

    % determine amplitude normalization factor
    normalizationFactor = sqrt(maximumSampleFrequency / sum(signal.^2));

    % normalize siganal to unity characteristic amplitude
    signal = normalizationFactor * signal;

    % handle odd length data
    if oddDataLength == 1,
      signal = signal(1 : end - 1);
    end

    % record injection properties
    properties.time = centralTime;
    properties.frequency = centralFrequency;
    properties.duration = 2 * sqrt(pi) * duration;
    properties.bandwidth = 1 / duration;
    properties.amplitude = 1;
    properties.phase = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                 inspiral                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'inspiral'

    % select random reference time [seconds]
    maximumSignalDuration = 4;
    minimumCentralTime = startTime + transientDuration + maximumSignalDuration;
    maximumCentralTime = stopTime - transientDuration - maximumSignalDuration;
    referenceTime = unifrnd(minimumCentralTime, maximumCentralTime);

    % select random component masses [solar mass]
    componentMass1 = unifrnd(1, 3);
    componentMass2 = unifrnd(1, 3);

    % select random reference phase [radians]
    referencePhase = unifrnd(0, 2 * pi);

    % normalize to horizon distance of unity [megaparsec]
    distance = 1;

    % reference frequency [Hz]
    referenceFrequency = 150;

    % fundamental constants
    speedOfLight = 2.99792458e10;
    gravitationalConstant = 6.6742e-8;

    % astrophysical constants
    solarMass = 1.98844e33;
    megaParsec = 3.0856775807e24;

    % convert to geometric units
    time = speedOfLight * time;
    componentMass1 = gravitationalConstant * solarMass * componentMass1 / ...
                     speedOfLight^2;
    componentMass2 = gravitationalConstant * solarMass * componentMass2 / ...
                     speedOfLight^2;
    referenceFrequency = referenceFrequency / speedOfLight;
    referenceTime = speedOfLight * referenceTime;
    referencePhase = referencePhase;
    distance = distance * megaParsec;

    % determine mass parameters
    totalMass = componentMass1 + componentMass2;
    reducedMass = componentMass1 * componentMass2 / totalMass;
    chirpMass = totalMass^(2/5) * reducedMass^(3/5);

    % time from reference to coalescence
    coalescenceTime = 5 * 256^(-1) * pi^(-8/3) * chirpMass^(-5/3) * ...
                      referenceFrequency^(-8/3);

    % dimensionless time parameter
    dimensionlessTime = 1 - (time - referenceTime) / coalescenceTime;

    % dimensionless time of last stable orbit
    dimensionlessEndTime = 6 * (pi * totalMass * referenceFrequency)^(8/3);

    % amplitude factor
    amplitude = 576^(1/4) * pi^(2/3) * chirpMass^(5/4) ./ ...
                (distance * (5 * coalescenceTime * dimensionlessTime).^(1/4));

    % accumulated phase
    accumulatedPhase = 16 * pi * referenceFrequency * coalescenceTime .* ...
                 (1 - dimensionlessTime.^(5/8)) / 5 + referencePhase;

    % generate signal
    signal = amplitude .* cos(accumulatedPhase);

    % truncate waveform at last stable orbit
    signal(find(dimensionlessTime < dimensionlessEndTime)) = 0;

    % convert back to standard units
    time = time / speedOfLight;
    referenceFrequency = referenceFrequency * speedOfLight;
    referenceTime = referenceTime / speedOfLight;
    referencePhase = referencePhase;
    coalescenceTime = coalescenceTime / speedOfLight;
    distance = distance / megaParsec;

    % record injection properties
    properties.time = referenceTime;
    properties.frequency = referenceFrequency;
    properties.duration = 2 * coalescenceTime;
    properties.bandwidth = 2 * referenceFrequency;
    properties.amplitude = distance;
    properties.phase = referencePhase;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                 ringdown                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'ringdown'

    % select random reference time [seconds]
    maximumSignalDuration = 1;
    minimumCentralTime = startTime + transientDuration + maximumSignalDuration;
    maximumCentralTime = stopTime - transientDuration - maximumSignalDuration;
    referenceTime = unifrnd(minimumCentralTime, maximumCentralTime);

    % select random central frequency [Hz]
    minimumCentralFrequency = tiling.minimumFrequency;
    maximumCentralFrequency = tiling.maximumFrequency;
    centralFrequency =  2.^(unifrnd(log2(minimumCentralFrequency), ...
                                    log2(maximumCentralFrequency)));

    % select random quality factor []
    maximumDuration = 1;
    minimumQualityFactor = 2;
    maximumQualityFactor = pi * maximumDuration * centralFrequency;
    qualityFactor =  2.^(unifrnd(log2(minimumQualityFactor), ...
                                 log2(maximumQualityFactor)));

    % select random phase [radians]
    minimumPhase = 0;
    maximumPhase = 2 * pi;
    phase = unifrnd(minimumPhase, maximumPhase);

    % normalize for unity characteristic amplitude
    amplitude = 1;

    % exponentially decaying envelope
    envelope = exp(-pi * centralFrequency * (time - referenceTime) / ...
                   qualityFactor);
    envelope(find(time < referenceTime)) = 0;

    % sinusoidal oscillation
    sinusoid = sin(2 * pi * centralFrequency * (time - referenceTime) + phase);

    % generate signal
    signal = amplitude * envelope .* sinusoid;

    % determine amplitude normalization factor
    normalizationFactor = sqrt(maximumSampleFrequency / sum(signal.^2));

    % normalize siganal to unity characteristic amplitude
    signal = normalizationFactor * signal;

    % record injection properties
    properties.time = referenceTime;
    properties.frequency = centralFrequency;
    properties.duration = qualityFactor / (pi * centralFrequency);
    properties.bandwidth = pi * centralFrequency / qualityFactor;
    properties.amplitude = 1;
    properties.phase = phase;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                 unknown                                    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  otherwise

    error(['ERROR: unknown waveform: ''' waveform '''']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         end switch on waveform type                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       estimate signal to noise ratios                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize injection data structure
data = cell(numberOfChannels, 1);

% initialize signal to noise estimates
signalToNoiseRatios = zeros(numberOfChannels, 1);

% begin loop over channels
for channelNumber = 1 : numberOfChannels,
  
  % extract sample frequency for channel
  sampleFrequency = sampleFrequencies(channelNumber);

  % length of noise data
  dataLength = length(noise{channelNumber});
  halfDataLength = dataLength / 2 + 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                       resample and reshape signal                          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % channel time vector
  channelTime = startTime + (0 : dataLength - 1) / sampleFrequency;
  
  % interpolate and reshape injection signal to noise sample frequency and shape
  data{channelNumber} = reshape(interp1(time, signal, channelTime, 'cubic'), ...
                                size(noise{channelNumber}));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         estimate noise spectrum                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % two-sided and one-sided fourier transform lengths
  transformLength = sampleFrequency / spectrumResolution;
  overlapLength = transformLength / 2;
  halfTransformLength = transformLength / 2 + 1;

  % window for spectral estimation
  window = hanning(transformLength);

  % list of fourier transform start indices
  startIndices = 1 : overlapLength : dataLength - transformLength + 1;
  numberOfTransforms = length(startIndices);

  % initialize noise spectrum estimate
  noiseSpectrum = zeros(transformLength, numberOfTransforms);

  % begin loop over fourier transforms
  for transformNumber = 1 : numberOfTransforms,

    % extract noise segment to transform
    startIndex = startIndices(transformNumber);
    stopIndex = startIndex + transformLength - 1;
    noiseSegment = noise{channelNumber}(startIndex : stopIndex);
    noiseSegment = window .* noiseSegment(:);
    noiseSpectrum(:, transformNumber) = fft(noiseSegment, transformLength);

  % end loop over fourier transforms
  end

  % extract one sided noise spectra
  noiseSpectrum = noiseSpectrum(1 : halfTransformLength, :);

  % calculate unbiased median estimate of noise spectrum
  noiseSpectrum = 2 * median(abs(noiseSpectrum).^2, 2) / ...
                  (log(2) * sampleFrequency * sum(window.^2));
  
  % interpolate one sided noise spectrum to half data length
  resolutionRatio = dataLength / transformLength;
  noiseSpectrum = ones(resolutionRatio, 1) * noiseSpectrum';
  noiseSpectrum = noiseSpectrum(:);
  noiseSpectrum = noiseSpectrum(resolutionRatio / 2 : ...
                                end - resolutionRatio / 2)';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      estimate signal to noise ratio                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % calculate signal spectrum
  signalSpectrum = abs(fft(data{channelNumber}(:)) / sampleFrequency).^2;
  signalSpectrum = signalSpectrum(1 : halfDataLength);

  % corresponding frequency vector
  signalFrequency = 0 : sampleFrequency / dataLength : sampleFrequency / 2;
  
  % bandlimit signal to noise ratio estimation
  inBandIndices = find((signalFrequency > min(spectrumBand)) & ...
                       (signalFrequency < max(spectrumBand)));
  signalSpectrum = signalSpectrum(inBandIndices);
  noiseSpectrum = noiseSpectrum(inBandIndices);

  % scale factor for desired signal to noise ratio
  snrIntegrand = 4 * signalSpectrum(:) ./ noiseSpectrum(:);
  frequencyResolution = sampleFrequency / dataLength;
  signalToNoiseRatios(channelNumber) = ...
      sqrt(sum(snrIntegrand * frequencyResolution));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     end estimate signal to noise ratios                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           inject signal into noise                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine injection scale factor
scaleFactor = signalToNoiseRatio / sqrt(sum(signalToNoiseRatios.^2));

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  % inject signal into detector noise
  data{channelNumber} = scaleFactor * data{channelNumber} + ...
                        noise{channelNumber};

% end loop over channels
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      construct injection description                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize result structures
injection = cell(numberOfChannels, 1);

% begin loop over channels
for channelNumber = 1 : numberOfChannels,
  
  % insert injection properties into injection description structure
  injection{channelNumber} = properties;
  injection{channelNumber}.snr = signalToNoiseRatio;
  switch lower(waveform)
    case 'inspiral'
      injection{channelNumber}.amplitude = ...
          injection{channelNumber}.amplitude / scaleFactor;
    otherwise
      injection{channelNumber}.amplitude = ...
          injection{channelNumber}.amplitude * scaleFactor;
   end

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                return results                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
