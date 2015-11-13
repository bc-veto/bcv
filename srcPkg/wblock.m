function varargout = ...
        wblock(blockStartTime, tiling, eventTime, ...
               parameters, frameCache, outputFiles, debugLevel)
% WBLOCK Omega Pipeline function to analyze a single block of data
%
% usage: wblock(blockStartTime, tiling, eventTime, ...
%               parameterFile, frameCacheFile, outputFiles, debugLevel);
%
%   blockStartTime        gps block start time
%   tiling                q tiling structure
%   eventTime             eventTime to focus search around ([])
%   parameters            parameter structure
%   frameCache            frame cache structure
%   outputFiles           structure of output file names
%   debugLevel            verboseness level of debug output
%
% If given output arguments, WBLOCK will output the following fields:
%
%   triggersThresholded   thresholded triggers
%   triggersDownselected  downselected triggers
%   triggersVetoed        vetoed triggers
%   triggersClustered     clustered triggers
%   triggersCoincident    coincident triggers
%   event                 bayesian event
%   skymap                bayesian skymap
%   transforms            transform structure
%
% The specified debugLevel controls the amount of detail in the output log.
% A debugLevel of unity is assumed by default.
%
% See also WPARAMETERS, WTILE, WREADDATA, WRESAMPLE, WCONDITION, WTRANSFORM,
% WTHRESHOLD, WSELECT, and WWRITEEVENTS.

% Authors:
% Jameson Rollins <jrollins@phys.columbia.edu>
% Antony Searle <antony.searle@anu.edu.au>
% Shourov K. Chatterji <shourov@ligo.caltech.edu>
% Leo C. Stein <lstein@ligo.mit.edu>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         parse command line arguments                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 7, nargin));

% apply default arguments
if (nargin < 3) || isempty(eventTime),
  eventTime = [];
end
if (nargin < 4) || isempty(parameters),
  parameters = wparameters;
end
if (nargin < 5) || isempty(frameCache),
  frameCache = loadframecache;
end
if (nargin < 6) || isempty(outputFiles),
  outputFiles = [];
end
if (nargin < 7) || isempty(debugLevel),
  debugLevel = 1;
end

% ensure numeric time and debug level
if ischar(eventTime),
  eventTime = str2double(eventTime);
end
if ischar(debugLevel),
  debugLevel = str2double(debugLevel);
end

% initialize output structures
varargout = cell(1,nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               read parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract parameters from parameters structure
channelNames = parameters.channelNames;
frameTypes = parameters.frameTypes;

analysisMode = parameters.analysisMode;
blockDuration = parameters.blockDuration;
sampleFrequency = parameters.sampleFrequency;
falseEventRate = parameters.falseEventRate;
timeShifts = parameters.timeShifts;
stateNames = parameters.stateNames;
stateTypes = parameters.stateTypes;
stateMasks = parameters.stateMasks;
injectionNames = parameters.injectionNames;
injectionTypes = parameters.injectionTypes;
injectionFactors = parameters.injectionFactors;
injectionTimeShifts = parameters.injectionTimeShifts;
doubleWhiten = parameters.doubleWhiten;
outlierFactor = parameters.outlierFactor;
maximumSignificants = parameters.maximumSignificants;
maximumTriggers = parameters.maximumTriggers;
durationInflation = parameters.durationInflation;
bandwidthInflation = parameters.bandwidthInflation;
coincidenceNumber = parameters.coincidenceNumber;
maximumCoincidents = parameters.maximumCoincidents;
triggerFields = parameters.triggerFields;
triggerFormat = parameters.triggerFormat;
numberOfChannels = parameters.numberOfChannels;
injectionChannels = parameters.injectionChannels;

% clustering parameters
applyClustering = parameters.applyClustering;
clusterMethod = parameters.clusterMethod;
clusterParameter1 = parameters.clusterParameter1;
clusterParameter2 = parameters.clusterParameter2;
clusterParameter3 = parameters.clusterParameter3;
distanceMetric = parameters.distanceMetric;
writeClusters = parameters.writeClusters;

% veto parameters
applyVeto = parameters.applyVeto;
falseVetoRate = parameters.falseVetoRate;
uncertaintyFactor = parameters.uncertaintyFactor;
correlationFactor = parameters.correlationFactor;
vetoDurationFactor = parameters.vetoDurationFactor;
vetoBandwidthFactor = parameters.vetoBandwidthFactor;
maximumConsistents = parameters.maximumConsistents;

% sky coordinate parameters
skyPosition = parameters.skyPosition;
skyCoordinateSystem = parameters.skyCoordinateSystem;

% variable to say if coincident events are found
coincidentsFound = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                determine times                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force times into alignment with data samples
blockStartTime = floor(blockStartTime * sampleFrequency) / sampleFrequency;

% stop time of block
blockStopTime = blockStartTime + blockDuration;

% force time shift alignment with sample interval
timeShifts = round(timeShifts * sampleFrequency) / sampleFrequency;

% report status
wlog(debugLevel, 1, '  block start time:        %#13.2f\n', blockStartTime);
wlog(debugLevel, 1, '  block stop time:         %#13.2f\n', blockStopTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           transform sky coordinate                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$%%%%%%%%%%%%%%%%%%%%%%

% if the skyPosition parameter is not empty, determine coordinate
if ~isempty(skyPosition),

  % gps center time of block
  blockCenterTime = blockStartTime + (blockDuration / 2);
  wlog(debugLevel, 2, '  block center time:       %#13.2f\n', blockCenterTime);

  % convert sky coordinate into geocentric coordinate
  coordinate = ...
      wconvertskycoordinates(skyPosition, blockCenterTime, ...
                             skyCoordinateSystem, 'geocentric');

  wlog(debugLevel, 1, '  sky coordinate:          [%.3f %.3f]\n', ...
       coordinate(1), coordinate(2));

else

    % otherwise output empty coordinate
    coordinate = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             read detector state                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if detector state channels are requested
if ~isempty(stateNames),

  % report status
  wlog(debugLevel, 1, '  reading detector state\n');

  % read detector data
  [stateData, stateSampleFrequencies] = ...
      wreaddata(frameCache, stateNames, stateTypes, blockStartTime, ...
                blockStopTime, timeShifts, debugLevel);

  % bool array for tracking valid channels
  validChannels = true(numberOfChannels,1);
  
  % loop over channels
  for channelNumber = 1 : numberOfChannels,

    % check for missing data
    if stateSampleFrequencies(channelNumber) == 0,
        wlog(debugLevel, 1, '    %-21s        state read error\n', ...
                [stateNames{channelNumber} ':']);
        validChannels(channelNumber) = false;
        continue;
    end

    % check for missing data
    if any(stateData{channelNumber} < 0),
        wlog(debugLevel, 1, '    %-21s        state data negative??\n', ...
                [stateNames{channelNumber} ':']);
        validChannels(channelNumber) = false;
        continue;
    end

    % check state
    if any(bitand(stateData{channelNumber}, stateMasks(channelNumber)) ~= ...
           stateMasks(channelNumber)),
        wlog(debugLevel, 1, '    %-21s        detector not in requested state\n', ...
                [stateNames{channelNumber} ':']);
        validChannels(channelNumber) = false;
        continue;
    end

  end

  % remove channels that didn't pass checks  
  channelNames = {channelNames{validChannels}}';
  frameTypes = {frameTypes{validChannels}}';
  timeShifts = (timeShifts(validChannels))';
  stateNames = {stateNames{validChannels}}';
  stateTypes = {stateTypes{validChannels}}';
  stateMasks = (stateMasks(validChannels))';
  injectionNames = {injectionNames{validChannels}}';
  injectionTypes = {injectionTypes{validChannels}}';

  % remove channels from output files list as well
  fieldNames = fieldnames(outputFiles);
  fieldNames = fieldNames(~strcmp(fieldNames, 'EVENTS') & ...
                          ~strcmp(fieldNames, 'SKYMAP_BASE'));
  for fieldNumber = 1 : length(fieldNames),
      outputFiles.(fieldNames{fieldNumber}) = ...
          {outputFiles.(fieldNames{fieldNumber}){validChannels}};
  end

  numberOfChannels = length(channelNames);

  % free detector state data memory
  clear stateData;
  clear stateSampleFrequencies;

% end test for requested detector state
end

% error if there is no data after state checks
if isempty(channelNames)
  wlog(debugLevel, 1, '    no valid data to analyze\n');
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             read detector data                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  reading detector data\n');

% read detector data
[data, sampleFrequencies] = ...
    wreaddata(frameCache, channelNames, frameTypes, blockStartTime, ...
              blockStopTime, timeShifts, debugLevel);

% test for missing detector data
if any(sampleFrequencies == 0),
  error('error reading detector data');
end

% standardize channel names
channelNames = strrep(channelNames, ';', ':');

% report significant tile numbers
for channelNumber = 1 : length(channelNames),
  wlog(debugLevel, 1, '    %-21s        loaded\n', ...
          [channelNames{channelNumber} ':']);
end

% determine network string
networkString = [];
for channelNumber = 1 : length(channelNames),
   networkString = sprintf('%s%s', networkString, ...
                           channelNames{channelNumber}(1:2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            read injection data                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if any injections are requested
if any(injectionChannels),

  % report status
  wlog(debugLevel, 1, '  reading injection data\n');

  % read injection data
  [injectionData, injectionSampleFrequencies] = ...
      wreaddata(frameCache, injectionNames, injectionTypes, blockStartTime, ...
                blockStopTime, injectionTimeShifts, debugLevel);

  % test for missing injection data
  if any(injectionSampleFrequencies(injectionChannels) == 0),
    error('error reading injection data');
  end

  % test for inconsistent injection data
  if any(injectionSampleFrequencies(injectionChannels) ~= ...
         sampleFrequencies(injectionChannels)),
    error('inconsistent detector and injection data');
  end

  % add injection data to detector data
  for channelNumber = injectionChannels,
    data{channelNumber} = data{channelNumber} + ...
        injectionFactors(channelNumber) * injectionData{channelNumber};
  end

  % free injection memory
  clear injectionData;
  clear injectionSampleFrequencies;

% end test for requested injections
end

% clear the frame cache, since all data is loaded
clear frameCache;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               resample data                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  resampling data\n');

% resample data
data = wresample(data, sampleFrequencies, sampleFrequency);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               condition data                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  conditioning data\n');

% condition data
[data, coefficients] = wcondition(data, tiling, doubleWhiten);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             transform and threshold                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup thresholds ranges
thresholdReferenceTime = [];
thresholdTimeRange = [];
thresholdFrequencyRange = [];
thresholdQRange = [];

% for coherent mode the transforming and thresholding have to be done
% separately, since the coherent threshold is based on the significance of
% reference channels (for 'signal' and 'null'), and therefore all transform
% channels have to be available before thresholding.
if strcmpi(analysisMode, 'coherent'),
  wlog(debugLevel, 1, '  transforming data\n');
  transforms = ...
      wtransform(data, tiling, outlierFactor, ...
                 analysisMode, channelNames, coefficients, coordinate);

  wlog(debugLevel, 1, '  thresholding data\n');
  triggersThresholded = ...
      wthreshold(transforms, tiling, blockStartTime, falseEventRate, ...
                 thresholdReferenceTime, thresholdTimeRange, ...
                 thresholdFrequencyRange, thresholdQRange, ...
                 maximumSignificants, analysisMode, falseVetoRate, ...
                 uncertaintyFactor, correlationFactor, debugLevel);

% if the mode is not coherent, then the transforming and thresholding can happen
% at once
else
  wlog(debugLevel, 1, '  transforming and thresholding data\n');
  triggersThresholded = ...
      wtransformandthreshold(data, blockStartTime, tiling, ...
               outlierFactor, falseEventRate, ...
               analysisMode, channelNames, coefficients, coordinate, ...
               thresholdReferenceTime, thresholdTimeRange, ...
               thresholdFrequencyRange, thresholdQRange, ...
               maximumSignificants, falseVetoRate, uncertaintyFactor, ...
               correlationFactor, debugLevel);

end

% set initial triggers
triggers = triggersThresholded;

% report significant tile numbers
for channelNumber = 1 : length(triggersThresholded),
  wlog(debugLevel, 1, '    %-21s %6u tiles\n', ...
          [triggersThresholded{channelNumber}.channelName ':'], ...
          length(triggersThresholded{channelNumber}.time));
end

% free data memory if we're not outputing it, and we're not in bayesian mode
% (bayesian search needs the raw data)
if nargout < 8 && ~strcmp(analysisMode, 'bayesian')
  clear data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        select significant triggers                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  selecting significant triggers\n');

% excluding overlapping significant tiles
triggersDownselected = ...
    wselect(triggers, durationInflation, bandwidthInflation, ...
               maximumTriggers, debugLevel);

% replace triggers
triggers = triggersDownselected;

% report downselected tile numbers
for channelNumber = 1 : length(triggersDownselected),
  wlog(debugLevel, 1, '    %-21s %6u tiles\n', ...
          [triggersDownselected{channelNumber}.channelName ':'], ...
          length(triggersDownselected{channelNumber}.time));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           write single detector triggers                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outputFiles) && isfield(outputFiles,'DOWNSELECT') && ...
        ~isempty(outputFiles.DOWNSELECT),

  % report status
  wlog(debugLevel, 1, '  writing selected triggers\n');

  % write triggers
  wwriteevents(triggersDownselected, outputFiles.DOWNSELECT, ...
               triggerFields, triggerFormat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          veto inconsistent triggers                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if veto is requested and analysis mode is coherent,
if applyVeto && strcmpi(analysisMode, 'coherent'),

  % report status
  wlog(debugLevel, 1, '  applying veto\n');

  % apply null stream veto
  triggersVetoed = ...
      wveto(triggers, durationInflation, bandwidthInflation, ...
            vetoDurationFactor, vetoBandwidthFactor, ...
            maximumConsistents, debugLevel);

  % replace triggers
  triggers = triggersVetoed;

  % report consistent tile numbers
  for channelNumber = 1 : length(triggersVetoed),
    wlog(debugLevel, 1, '    %-21s %6u tiles\n', ...
            [triggersVetoed{channelNumber}.channelName ':'], ...
            length(triggersVetoed{channelNumber}.time));
  end
  
% else make empty veto trigger structure
else
    
  triggersVetoed = [];
    
% end test for apply veto
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                cluster triggers                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if clustering is requested,
if applyClustering || strcmp(analysisMode, 'bayesian'),

  % report status
  wlog(debugLevel, 1, '  applying clustering\n');

  % apply clusterings
  [triggersClustered, clusters] = ...
  wcluster(triggers, clusterMethod, ...
           clusterParameter1, clusterParameter2, clusterParameter3, ...
           distanceMetric, durationInflation, bandwidthInflation, ...
           debugLevel);

  % replace triggers
  triggers = triggersClustered;

  % report cluster numbers
  for channelNumber = 1 : length(triggersClustered),
    wlog(debugLevel, 1, '    %-21s %6u clusters\n', ...
            [triggersClustered{channelNumber}.channelName ':'], ...
            length(clusters{channelNumber}.time));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         write cluster triggers                             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~isempty(outputFiles) && isfield(outputFiles,'CLUSTER') && ...
        ~isempty(outputFiles.CLUSTER),

    % report status
    wlog(debugLevel, 1, '  writing cluster triggers\n');

    % write triggers or clusters
    if writeClusters,
      wwriteevents(clusters, outputFiles.CLUSTER, ...
                   triggerFields, triggerFormat);
    else
      wwriteevents(triggersClustered, outputFiles.CLUSTER, ...
                   triggerFields, triggerFormat);
    end

  end

% else make empty cluster trigger structure
else

  triggersClustered = [];

% end test for apply clustering
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              coincide triggers                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply coincidence if requested
if (coincidenceNumber > 1) && (numberOfChannels > 1),

  % report status
  wlog(debugLevel, 1, '  coinciding triggers\n');

  % apply coincidence
  triggersCoincident = ...
      wcoincide(triggers, coincidenceNumber, durationInflation, ...
                bandwidthInflation, maximumCoincidents, debugLevel);

  % replace triggers
  triggers = triggersCoincident;

  % cull coincident clusters based on coincident triggers
  if exist('clusters','var'),
    clustersCoincident = ...
        wcoincidecluster(triggersCoincident, clusters, debugLevel);

    % replace clusters
    clusters = clustersCoincident;
  end

  % note if any coincidents were found
  if ~isempty(triggersCoincident{1}.amplitude),
    coincidentsFound = true;
  end

  % report coincidence numbers
  for channelNumber = 1 : length(triggersCoincident),
    wlog(debugLevel, 1, '    %-21s %6u coincidents\n', ...
            [triggersCoincident{channelNumber}.channelName ':'], ...
            length(triggersCoincident{channelNumber}.time));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         write coincident triggers                          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~isempty(outputFiles) && isfield(outputFiles,'COINCIDE') && ...
          ~isempty(outputFiles.COINCIDE),

    % report status
    wlog(debugLevel, 1, '  writing coincident triggers\n');

    % write triggers or clusters
    if applyClustering && writeClusters,
      wwriteevents(clustersCoincident, outputFiles.COINCIDE, ...
                   triggerFields, triggerFormat);
    else
      wwriteevents(triggersCoincident, outputFiles.COINCIDE, ...
                   triggerFields, triggerFormat);
    end

  end

% else make empty coincident trigger structure
else

  triggersCoincident = [];
  clustersCoincident = [];

% end test for apply coincidence
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              bayesian followup                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply followup if requested and coincident tiles and clusters have been found
if strcmp(analysisMode, 'bayesian') && coincidentsFound,

  % report status
  wlog(debugLevel, 1, '  performing followup on significant event\n');

  % find candidate cluster
  candidateCluster = ...
      wfindcandidate(clustersCoincident, triggersClustered, ...
                     eventTime, debugLevel);

  % apply followup to find most significant event and event skymap
  [event, skymap] = ...
      wfollowup(data, coefficients, tiling, blockStartTime, coordinate, ...
                candidateCluster, channelNames, parameters, debugLevel);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                       write additional event fields                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % write in network string
  for channelNumber = 1 : numberOfChannels,
    if channelNumber == 1,
      event.network =  channelNames{channelNumber}(1:2);
    else
      event.network = sprintf('%s,%s',event.network,channelNames{channelNumber}(1:2));
    end
  end

  % write block times
  event.blockStartTime = blockStartTime;
  event.blockStopTime = blockStopTime;
  % write total livetime
  event.livetime = blockDuration - 2 * tiling.transientDuration;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          write events and skymap                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % write out the event and skymap, if they are not empty (ie. pass threshold)
  if ~isempty(event) && ~isempty(outputFiles) ...
          && isfield(outputFiles,'EVENTS') && ...
          ~isempty(outputFiles.EVENTS),

    % report status
    wlog(debugLevel, 1, '  writing event\n');

    % event fields
    eventFields = {'network', ...
                   'blockStartTime', 'blockStopTime', 'livetime', ...
                   'time', 'frequency', ...
                   'duration', 'bandwidth', ...
                   'logSignal', 'logGlitch', ...
                   'modeTheta', 'modePhi', ...
                   'probSignal', 'probGlitch'};

    % write events
    wwriteevents(event, outputFiles.EVENTS, eventFields, triggerFormat)

  end

  % write out the skymap
  if ~isempty(skymap) && ~isempty(outputFiles) ...
          && isfield(outputFiles,'SKYMAP_BASE') && ...
          ~isempty(outputFiles.SKYMAP_BASE),

    % report
    wlog(debugLevel, 1, '  writing skymap\n');

    % write skymap for event time
    wwriteskymap(skymap, ...
                 strrep(outputFiles.SKYMAP_BASE, ...
                        '@TIME@', sprintf('%#020.9f', event.time)));

  end

% else make empty event and skymap
else

  event = [];
  skymap = [];

% end test for apply followup
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         process output arguments                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 1,
    varargout{1} = channelNames;
end
if nargout > 1,
    varargout{1} = triggersThresholded;
    varargout{2} = triggersDownselected;
    varargout{3} = triggersVetoed;
    varargout{4} = triggersClustered;
    varargout{5} = triggersCoincident;
    varargout{6} = event;
    varargout{7} = skymap;
end
if nargout >= 8,
    % output the transform of the data, for plotting spectrograms
    if ~exist('transforms','var'),
      transforms = ...
          wtransform(data, tiling, outlierFactor, analysisMode, ...
                     channelNames, coefficients, coordinate);
    end
    varargout{8} = transforms;
end
