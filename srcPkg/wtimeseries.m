function handles = wtimeseries(data, tiling, startTime, referenceTime, ...
                               timeRange, channelNames)
% WTIMESERIES Display time series data from Q transform analysis
%
% WTIMESERIES displays raw, high pass filtered, or whitened time series
% data used in Q transform analysis.  A separate figure is produced for
% each channel under analysis.
%
% usage:
%
%   handles = wtimeseries(data, tiling, startTime, referenceTime, ...
%                         timeRange, channelNames);
%
%     data             cell array of time domain data
%     tiling           q transform tiling structure
%     startTime        start time of data
%     referenceTime    reference time of plot
%     timeRange        vector range of relative times to plot
%     channelNames     cell array of channel names
%
%     handles          matrix of axis handles for each plot
%
% To determine the range of times to plot, WSPECTROGRAM requires the start time
% of the data under analysis, a reference time, and relative time range.  The
% relative time range should be specified as a two component vector, consisting
% of the minimum and maximum relative times to plot.  The default choice of all
% available data can be obtained by passing the empty matrix [].  Both the start
% time and reference time should be specified as absolute quantities, but the
% range of times to plot should be specified relative to the requested reference
% time.  The specified reference time is used as the time origin in the
% resulting spectrograms and is also reported in the title of each plot.
%
% WTIMESERIES expects time domain data.  To display the time series of
% conditioned data returned by WCONDITION or WSCANCONDITION, first use WIFFT to
% convert the one-sided frequency domain data to the time domain.
%
% If only one channel is requested, it is plotted in the current figure window.
% If more than one channel are requested, they are plotted in separate figure
% windows starting with figure 1.
%
% The optional cell array of channel names are used to label the resulting
% figures.
%
% WTIMESERIES returns a vector of axis handles with one handle per channel.
%
% See also WFFT, WIFFT, WCONDITION, WSCANCONDITION, WSPECTROGRAM, WEVENTGRAM,
% and WSCAN.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu

% $Id: wtimeseries.m 986 2008-08-14 21:04:56Z lstein $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            hard coded parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image boundary
imageLeft = 0.14;
imageWidth = 0.80;
imageBottom = 0.12;
imageHeight = 0.78;
imagePosition = [imageLeft imageBottom imageWidth imageHeight];

% time scales for labelling
millisecondThreshold = 0.5;
secondThreshold = 3 * 60;
minuteThreshold = 3 * 60 * 60;
hourThreshold = 3 * 24 * 60 * 60;
dayThreshold = 365.25 * 24 * 60 * 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arugments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 6, nargin));

% apply default arguments
if (nargin < 3) || isempty(startTime),
  startTime = 0;
end
if (nargin < 4) || isempty(referenceTime),
  referenceTime = 0;
end
if (nargin < 5) || isempty(timeRange),
  timeRange = [-Inf Inf];
end
if (nargin < 6) || isempty(channelNames),
  channelNames = [];
end

% force cell arrays
data = wmat2cell(data);
channelNames = wmat2cell(channelNames, ~isempty(channelNames));

% force one dimensional cell arrays
data = data(:);
channelNames = channelNames(:);

% determine number of channels
numberOfChannels = length(data);

% provide default channel names
if isempty(channelNames),
  channelNames = cell(numberOfChannels, 1);
  for channelNumber = 1 : numberOfChannels,
    channelNames{channelNumber} = ['Channel ' int2str(channelNumber)];
  end
end

% force time range to be monotonically increasing column vector
timeRange = unique(timeRange(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate channel names
if ~isempty(channelNames) && (length(channelNames) ~= numberOfChannels),
  error('channel names is inconsistent with number of transform channels');
end

% check for valid discrete Q-transform tile structure
if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
  error('The first argument is not a valid Q-transform tiling.');
end

% check for two component time range vector
if length(timeRange) ~= 2,
  error('Time range must be two component vector [tmin tmax].');
end

% validate data lengths
dataLength = tiling.duration * tiling.sampleFrequency;
for channelNumber = 1 : numberOfChannels,
  if length(data{channelNumber}) ~= dataLength,
    error('Data length is inconsistent with tiling');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          identify times to display                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relative time vector
relativeTime = (startTime - referenceTime) + ...
    (0 : (dataLength - 1)) / tiling.sampleFrequency;

% default start time for display is start time of available data
if timeRange(1) == -Inf,
  timeRange(1) = startTime - referenceTime;
end

% default stop time for display is stop time of available data
if timeRange(2) == +Inf,
  timeRange(2) = startTime - referenceTime + tiling.duration;
end

% validate requested time range
if (timeRange(1) < startTime - referenceTime) || ...
   (timeRange(2) > startTime - referenceTime + tiling.duration),
  error('requested time range exceeds available data');
end

% indices of times to display
plotIndices = find((relativeTime >= min(timeRange)) & ...
                   (relativeTime <= max(timeRange)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize handle vector
handles = zeros(numberOfChannels);

% loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              plot time series                              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % if plotting more than one figure
  if numberOfChannels > 1,

    % select figure to plot in
    figure(channelNumber);

  % continue
  end

  % reset figure
  clf;
  set(gca, 'FontSize', 16);

  % determine order of magnitude ammplitude scale
  amplitudeScale = 10^floor(log10(max(abs(data{channelNumber}(plotIndices)))));

  % plot time series
  if abs(diff(timeRange)) < millisecondThreshold,
    plot(relativeTime(plotIndices) * 1e3, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  elseif abs(diff(timeRange)) < secondThreshold,
    plot(relativeTime(plotIndices) * 1, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  elseif abs(diff(timeRange)) < minuteThreshold,
    plot(relativeTime(plotIndices) / 60, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  elseif abs(diff(timeRange)) < hourThreshold,
    plot(relativeTime(plotIndices) / 3600, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  elseif abs(diff(timeRange)) < dayThreshold,
    plot(relativeTime(plotIndices) / 86400, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  else
    plot(relativeTime(plotIndices) / 31557600, ...
         data{channelNumber}(plotIndices) / amplitudeScale, 'b-');
  end

  % set axis position
  set(gca, 'Position', imagePosition);

  % set paper position
  set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]);

  % disable coordinate grid
  grid off;

  % set axis range
  if abs(diff(timeRange)) < millisecondThreshold,
    set(gca, 'XLim', timeRange * 1e3);
  elseif abs(diff(timeRange)) < secondThreshold,
    set(gca, 'XLim', timeRange * 1);
  elseif abs(diff(timeRange)) < minuteThreshold,
    set(gca, 'XLim', timeRange / 60);
  elseif abs(diff(timeRange)) < hourThreshold,
    set(gca, 'XLim', timeRange / 3600);
  elseif abs(diff(timeRange)) < dayThreshold,
    set(gca, 'XLim', timeRange / 86400);
  else
    set(gca, 'XLim', timeRange / 31557600);
  end

  % set x axis properties
  if abs(diff(timeRange)) < millisecondThreshold,
    xlabel('Time [milliseconds]');
  elseif abs(diff(timeRange)) < secondThreshold,
    xlabel('Time [seconds]');
  elseif abs(diff(timeRange)) < minuteThreshold,
    xlabel('Time [minutes]');
  elseif abs(diff(timeRange)) < hourThreshold,
    xlabel('Time [hours]');
  elseif abs(diff(timeRange)) < dayThreshold,
    xlabel('Time [days]');
  else
    xlabel('Time [years]');
  end

  % set y axis properties
  ylabel(['Amplitude [10^{' int2str(log10(amplitudeScale)) '}]']);

  % set title properties
  titleString = sprintf('%s at %.3f', ...
                        channelNames{channelNumber}, referenceTime);
  titleString = strrep(titleString, '_', '\_');
  title(titleString);

  % set figure background color
  set(gca, 'Color', [1 1 1]);
  set(gcf, 'Color', [1 1 1]);
  set(gcf, 'InvertHardCopy', 'off');

  % append current axis handle to list of handles
  handles(channelNumber) = gca;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return to calling function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
