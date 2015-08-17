function coincidents = ...
    wcoincide(events, coincidenceCount, durationInflation, ...
              bandwidthInflation, maximumCoincidents, debugLevel)
% WCOINCIDE Keep significant events that are coincident in several
% detectors
%
% WCOINCIDE discards significant signal events unless they overlap in time 
% and frequency with significant events in all other channels, allowing for
% arrival time.
%
% WCOINCIDE returns a cell array of coincident significant event properties
%
% usage: coincidents = ...
%            wcoincide(events, coincidenceCount, durationInflation, ...
%                bandwidthInflation, maximumCoincidents, debugLevel);
%
%   events               cell array of event properties
%   coincidenceCount     number of channels an event must be consistent 
%                        in; if empty all channels must be conincident
%   durationInflation    multiplicative scale factor for veto duration
%   bandwidthInflation   multiplicative scale factor for veto bandwidth
%   maximumCoincidents   maximum allowable number of coincident events
%   debugLevel           verboseness of debug output
%
%   coincidents          cell array of coincident event properties
%
% WCOINCIDE expects at least two input channels, with events sorted in 
% order of decreasing significance.  WCOINCIDE returns the same number of 
% channels.
%
% The optional durationInflation and bandwidthInflation arguments are
% multiplicative scale factors that are applied to the duration and bandwidth of
% significant events prior to testing for overlap.  If not specified, these
% parameters both default to unity, such that the resulting events have unity
% time-frequency area.
%
% The optional maximumCoincidents argument provides a safety mechanism to limit
% the total number of consistent events returned by WCOINCIDE.  If this maximum
% number of coincident events is exceeded, an overflow flag is set, only the
% maximumCoincidents most significant coincident events are returned, and a
% warning is issued if debugLevel is set to unity or higher.  By default,
% maximumCoincidents is set to infinity and debugLevel is set to unity.
%
% WCOINCIDE both expects and returns a cell array of Q transform event structures
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
% See also WTILE, WCONDITION, WTRANSFORM, WTHRESHOLD, WSELECT, WVETO, 
% WEXAMPLE, and WSEARCH.

% Shourov K. Chatterji <shourov@ligo.caltech.edu>
% Antony C. Searle <acsearle@ligo.caltech.edu>

% $Id: wcoincide.m 2304 2009-09-04 20:44:35Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(5, 6, nargin));

if isempty(durationInflation),
  durationInflation = 1.0;
end
if isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
end
if isempty(maximumCoincidents),
  maximumCoincidents = Inf;
end
if (nargin < 4) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
events = wmat2cell(events);

% force one dimensional cell arrays
events = events(:);

% determine number of channels
numberOfChannels = length(events);

if isempty(coincidenceCount)
  coincidenceCount = numberOfChannels;
end

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
  error('at least two channels required for coincidence');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    initialize consistent event structure                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of significant event structures
coincidents = cell([numberOfChannels 1]);

for channelNumber = 1:numberOfChannels,

    % insert structure identification string
    coincidents{channelNumber}.id = 'Discrete Q-transform event structure';

    % propogate overflow flag
    coincidents{channelNumber}.overflowFlag = events{channelNumber}.overflowFlag;

    % output channel name
    channelNames{channelNumber} = ...
        regexprep(events{channelNumber}.channelName, ':.*$', ':COINCIDENT');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         get allowed intersite delays                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delays = wbaseline(channelNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           enforce coincidence                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate output
coincidentIndices = cell([numberOfChannels 1]);

% loop over channel to be culled
for outerChannelNumber = 1:numberOfChannels,
    % record the number of events in the channel to be culled
    numberOfEvents = length(events{outerChannelNumber}.time);
    % for each event in the channel to be culled
    for eventNumber = 1:numberOfEvents,
        % track if the event is still viable
        count = 1;
        % for every channel
        for innerChannelNumber = 1:numberOfChannels
            % if the event is still viable and the channel is not its
            % channel of origin
            if (count < coincidenceCount) && (innerChannelNumber ~= outerChannelNumber)
                % check that the event's centre lies within the
                % time-frequency rectangle of at least one event in the
                % other channel, allowing for the light travel time
                eventViable = ...
                    any(( ...
                        abs(events{outerChannelNumber}.frequency(eventNumber) - ...
                            events{innerChannelNumber}.frequency) ...
                        < ...
                        max(events{outerChannelNumber}.bandwidth(eventNumber), ...
                            events{innerChannelNumber}.bandwidth) * ...
                            bandwidthInflation / 2 ...
                    ) & ( ...
                        abs(events{outerChannelNumber}.time(eventNumber) - ...
                            events{innerChannelNumber}.time) ...
                        < ...
                        delays(outerChannelNumber, innerChannelNumber) + ...
                        max(events{outerChannelNumber}.duration(eventNumber), ...
                            events{innerChannelNumber}.duration) * ...
                            durationInflation / 2 ...
                    ));
		if eventViable
		  count = count + 1;
		end
            end
        end
        % if the event has passed all checks, declare it consistent
        if count >= coincidenceCount                
            coincidentIndices{outerChannelNumber} = [coincidentIndices{outerChannelNumber} eventNumber];              
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   prune the events to speed up the next round of checks   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract coincident signal events
    coincidents{outerChannelNumber} = ...
        wcopyevents(events{outerChannelNumber}, ...
                    coincidentIndices{outerChannelNumber});
    % new channel names
    coincidents{outerChannelNumber}.channelName = ...
        channelNames{outerChannelNumber};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               check for excessive number of consistent events                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for channelNumber = 1:numberOfChannels
    % determine number of coincident events
    numberOfCoincidents = length(coincidents{channelNumber}.time);

    % if maximum allowable number of coincident events is exceeded
    if numberOfCoincidents > maximumCoincidents,

        % issue warning
        wlog(debugLevel, 1, 'WARNING: %s: maximum number of coincidents exceeded.\n', ...
             coincidents{channelNumber}.channelName);

        % set overflow flag
        coincidents{channelNumber}.overflowFlag = 1;

        % extract most significant signal events
        coincidents{channelNumber} = ...
            wcopyevents(coincidents{channelNumber}, (1:maximumCoincidents));

    % otherwise continue
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           return consistent events                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
