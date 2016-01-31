function consistents = ...
    wveto(events, durationInflation, bandwidthInflation, ...
          vetoDurationFactor, vetoBandwidthFactor, ...
          maximumConsistents, debugLevel)
% WVETO Veto significant signal events that overlap significant null events
%
% WVETO vetoes significant signal events if they overlap in time and frequency
% with significant null events, where the vetoed time frequency area potentially
% depends on the significance of the null event.
%
% WVETO returns a cell array of vetoed significant event properties
%
% usage: consistents = ...
%            wveto(events, durationInflation, bandwidthInflation, ...
%                  vetoDurationFactor, vetoBandwidthFactor, ...
%                  maximumConsistents, debugLevel);
%
%   events               cell array of signal and null event properties
%   durationInflation    multiplicative scale factor for veto duration
%   bandwidthInflation   multiplicative scale factor for veto bandwidth
%   vetoDurationFactor   veto duration inflation scaling factor
%   vetoBandwidthFactor  veto bandwidth inflation scaling factor
%   maximumConsistents   maximum allowable number of consistent events
%   debugLevel           verboseness of debug output
%
%   consistents          cell array of consistent event properties
%
% WVETO expects two input channels.  The first channel is the signal channel,
% and the second is the null channel.  WVETO returns a single channel containing
% consistent signal events.
%
% The optional durationInflation and bandwidthInflation arguments are
% multiplicative scale factors that are applied to the duration and bandwidth of
% significant events prior to testing for overlap.  If not specified, these
% parameters both default to unity, such that the resulting events have unity
% time-frequency area.  The optional vetoDurationFactor and vetoBandwidthFactor
% arguments are scaling factors that define how the durationInflation and
% bandwidthInflation factors increase with the normalized energy of the veto
% event.  If not specified, these parameters both default to zero, such that the
% vetoed time frequency area does not depend on null event significance.
%
% The optional maximumConsistents argument provides a safety mechanism to limit
% the total number of consistent events returned by WVETO.  If this maximum
% number of consistent events is exceeded, an overflow flag is set, only the
% maximumConsistents most significant consistent events are returned, and a
% warning is issued if debugLevel is set to unity or higher.  By default,
% maximumConsistents is set to infinity and debugLevel is set to unity.
%
% WVETO assumes that the input signal and null event lists have already been
% sorted in order of decreasing significance, as returned by the WSELECT
% function.
%
% WVETO both expects and returns a cell array of Q transform event structures
% with one cell per channel.  The event structures contain the following
% required fields, which describe the properties of statistically significant
% events.  Additional fields such as amplitude, are optional and are retained
% along with the required fields.
%
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   duration             duration of tile [seconds]
%   bandwidth            bandwidth of tile [Hz]
%   normalizedEnergy     normalized energy of tile []
%   incoherentEnergy     reference normalized energy of tile []
%
% The event structures also contain the following flag which indicates if the
% maximum number of significant tiles, significant events, or consistent events
% was exceeded.
%
%   overflowFlag         boolean overflow flag
%
% See also WTILE, WCONDITION, WTRANSFORM, WTHRESHOLD, WSELECT, WEXAMPLE, and
% WSEARCH.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu

% $Id: wveto.m 2304 2009-09-04 20:44:35Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 7, nargin));

% apply default arguments
if (nargin < 2) || isempty(durationInflation),
  durationInflation = 1.0;
end
if (nargin < 3) || isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
end
if (nargin < 4) || isempty(vetoDurationFactor),
  vetoDurationFactor = 0.0;
end
if (nargin < 5) || isempty(vetoBandwidthFactor),
  vetoBandwidthFactor = 0.0;
end
if (nargin < 6) || isempty(maximumConsistents),
  maximumConsistents = Inf;
end
if (nargin < 7) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
events = wmat2cell(events);

% force one dimensional cell arrays
events = events(:);

% determine number of channels
numberOfChannels = length(events);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(events{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

% validate number of channels
if numberOfChannels < 2,
  wlog(debugLevel, 2, '  skipping veto for < 2 channels')
  consistents = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    initialize consistent event structure                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% signal and null channel numbers
signalChannel = 1;
nullChannel = 2;

% create empty cell array of significant event structures
consistents = cell(1);

% insert structure identification string
consistents{signalChannel}.id = 'Discrete Q-transform event structure';

% propogate overflow flag
consistents{signalChannel}.overflowFlag = events{signalChannel}.overflowFlag;

% output channel name
consistents{signalChannel}.channelName = ...
    regexprep(events{signalChannel}.channelName, ':.*$', ':CONSISTENT');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         find signal event boundaries                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vector of significant tile start times
minimumTimes{signalChannel} = events{signalChannel}.time - ...
    durationInflation * ...
    events{signalChannel}.duration / 2;

% vector of significant tile stop times
maximumTimes{signalChannel} = events{signalChannel}.time + ...
    durationInflation * ...
    events{signalChannel}.duration / 2;

% vector of significant tile lower frequencies
minimumFrequencies{signalChannel} = events{signalChannel}.frequency - ...
    bandwidthInflation * ...
    events{signalChannel}.bandwidth / 2;

% vector of significant tile upper frequencies
maximumFrequencies{signalChannel} = events{signalChannel}.frequency + ...
    bandwidthInflation * ...
    events{signalChannel}.bandwidth / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          find null event boundaries                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vector of significant tile start times
minimumTimes{nullChannel} = events{nullChannel}.time - ...
    max(durationInflation, vetoDurationFactor * ...
        sqrt(2 * events{nullChannel}.normalizedEnergy)) .* ...
    events{nullChannel}.duration / 2;

% vector of significant tile stop times
maximumTimes{nullChannel} = events{nullChannel}.time + ...
    max(durationInflation, vetoDurationFactor * ...
        sqrt(2 * events{nullChannel}.normalizedEnergy)) .* ...
    events{nullChannel}.duration / 2;

% vector of significant tile lower frequencies
minimumFrequencies{nullChannel} = events{nullChannel}.frequency - ...
    max(bandwidthInflation, vetoBandwidthFactor * ...
        sqrt(2 * events{nullChannel}.normalizedEnergy)) .* ...
    events{nullChannel}.bandwidth / 2;

% vector of significant tile upper frequencies
maximumFrequencies{nullChannel} = events{nullChannel}.frequency + ...
    max(bandwidthInflation, vetoBandwidthFactor * ...
        sqrt(2 * events{nullChannel}.normalizedEnergy)) .* ...
    events{nullChannel}.bandwidth / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  apply veto                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of signal events
numberOfSignals = length(minimumTimes{signalChannel});

% number of null events
numberOfNulls = length(minimumTimes{nullChannel});

% initialize list of consistent signal events
consistentIndices = 1 : numberOfSignals;

% begin loop over null events
for nullIndex = 1 : numberOfNulls,
  
  % find non-vetoed signal events
  keepIndices = ...
      find((minimumTimes{nullChannel}(nullIndex) > ...
            maximumTimes{signalChannel}(consistentIndices)) | ...
           (maximumTimes{nullChannel}(nullIndex) < ...
            minimumTimes{signalChannel}(consistentIndices)) | ...
           (minimumFrequencies{nullChannel}(nullIndex) > ...
            maximumFrequencies{signalChannel}(consistentIndices)) | ...
           (maximumFrequencies{nullChannel}(nullIndex) < ...
            minimumFrequencies{signalChannel}(consistentIndices)));
  
  % number of vetoed signal events
  numberOfVetoedSignals = length(consistentIndices) - length(keepIndices);
  
  % remove vetoed signal events
  consistentIndices = consistentIndices(keepIndices);

  % if no signals were vetoed and there are more null events than signals events
  if ((numberOfVetoedSignals == 0) && ...
      (numberOfNulls - nullIndex > length(consistentIndices))),

    % end loop over null events
    break
  
  % otherwise, continue
  end
  
% end loop over null events
end

% remove tested null events
if ((numberOfNulls > 0) && (nullIndex + 1 <= numberOfNulls)),
  minimumTimes{nullChannel} = ...
      minimumTimes{nullChannel}(nullIndex + 1 : end);
  maximumTimes{nullChannel} = ...
      maximumTimes{nullChannel}(nullIndex + 1 : end);
  minimumFrequencies{nullChannel} = ...
      minimumFrequencies{nullChannel}(nullIndex + 1 : end);
  maximumFrequencies{nullChannel} = ...
      maximumFrequencies{nullChannel}(nullIndex + 1 : end);
end

% number of remianing signal events
numberOfSignals = length(consistentIndices);

% intialize list of consistent events to keep
keepIndices = [];

% begin loop over remaining signal events
for consistentIndex = 1 : numberOfSignals,
  
  signalIndex = consistentIndices(consistentIndex);
  
  % determine if signal event overlaps any null events
  overlap = ...
      any((minimumTimes{signalChannel}(signalIndex) < ...
           maximumTimes{nullChannel}) & ...
          (maximumTimes{signalChannel}(signalIndex) > ...
           minimumTimes{nullChannel}) & ...
          (minimumFrequencies{signalChannel}(signalIndex) < ...
           maximumFrequencies{nullChannel}) & ...
          (maximumFrequencies{signalChannel}(signalIndex) > ...
           minimumFrequencies{nullChannel}));

  % if signal event does not overlap with any null event,
  if ~overlap,

    % append it to the list of consistent events to keep
    keepIndices = [keepIndices consistentIndex];

  % otherwise, continue
  end

% end loop over remaining signal events
end
  
% remove vetoed signal events
consistentIndices = consistentIndices(keepIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               check for excessive number of consistent events                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of consistent events
numberOfConsistents = length(consistentIndices);

% if maximum allowable number of consistent events is exceeded
if numberOfConsistents > maximumConsistents,

  % issue warning
  wlog(debugLevel, 1, 'WARNING: %s: maximum number of consistents exceeded.\n', ...
       consistents{signalChannel}.channelName);

  % set overflow flag
  consistents{signalChannel}.overflowFlag = 1;

  % indices of most significant consistent events
  consistentIndices = consistentIndices(1 : maximumConsistents);

% otherwise continue
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  construct list fo consistent signal events                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract consistent signal events
consistents{signalChannel} = ...
    wcopyevents(events{signalChannel}, consistentIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           return consistent events                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
