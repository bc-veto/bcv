function parameters = wparameters(parameterFile, debugLevel)
% WPARAMETERS Read WSEARCH parameter file
%
% WPARAMETERS reads the parameters from the specifiec parameter file for
% use by the WSEARCH, WEVENT, and WPROPERTIES functions.
%
% usage: parameters = qarameters(parameterFile, debugLevel)
%
%   parameterFile       path name of parameter file to read
%   debugLevel          verboseness level of debug output
%
%   parameters          parameter structure
%
% By default, WPARAMETERS assumes the parameterFile name './parameters.txt'
% and a debugLevel of 1.
%
% The parameter file syntax conists of a field name and value separate by a
% colon.  The value should be an expression evaluable by Matlab, and can
% include scalars, vectors, matrices, or cell arrays.
%
%   fieldName:             fieldValue
%
% Comments are delineated by '#', '%', or '//', and extend to the end of
% the current line.  Blank lines are allowed, and are simply ignored.
%
% The following fields are required:
%
%   channelNames:
%   frameTypes:
%
% The following fields are optional and have the listed default value:
%
%   analysisMode:          'independent'
%   sampleFrequency:       4096
%   qRange:                [sqrt(11) 100]
%   frequencyRange:        [48 Inf]
%   maximumMismatch:       0.2
%   falseEventRate:        1
%   blockDuration:         64 (seconds)
%   conditionDuration:     (value of blockDuration)
%   timeShifts:            0
%   injectionNames:        'NONE'
%   injectionTypes:        'NONE'
%   injectionFactors:      0
%   injectionTimeShifts:   0
%   highPassCutoff:        (determined from tiling)
%   lowPassCutoff:         (determined from tiling)
%   whiteningDuration:     (determined from tiling)
%   transientFactor:       4
%   doubleWhiten:          1
%   extraBlockOverlap:     0
%   outlierFactor:         2.0
%   maximumSignificants:   1e5
%   maximumTriggers:       1e3
%   durationInflation:     1.0
%   bandwidthInflation:    1.0
%   coincidenceNumber:     0
%   maximumCoincidents:    Inf
%   triggerFields:         (depends on analysis mode and clustering)
%   triggerFormat:         'txt'
%   randomSeed:            sum(1e6 * clock)
%
% The following optional fields are used to read detector state data.
%
%   stateNames:            
%   stateTypes:            
%   stateMasks:            
%
% The following optional fields related to clustering.
%
%   applyClustering:       0
%   clusterMethod:         'density'
%   clusterRadius:         4.0                   (for density clustering)
%   clusterDensity:        3.0                   (for density clustering)
%   clusterSingles:        1                     (for density clustering)
%   clusterLinkage:        'single'              (for hierarchical clustering)
%   clusterCriterion:      'distance'            (for hierarchical clustering)
%   clusterThreshold:      4.0                   (for hierarchical clustering)
%   distanceMetric:        'integratedMismatch'
%   writeClusters:         0
%
% The following optional fields are for targeted searches
%
%   skyPosition:           []
%   skyCoordinateSystem:   'equatorial'
%
% The following fields apply only to 'coherent' anlaysis mode.
%
%   applyVeto:             1
%   falseVetoRate:         0.0
%   uncertaintyFactor:     0.0
%   correlationFactor:     0.0
%   vetoDurationFactor:    0.5
%   vetoBandwidthFactor:   0.5
%   maximumConsistents:    1e3
%
% The following fields apply only in 'bayesian' followup analysis mode.
%
%   maxFollowTriggers      5
%   xCoherentCheck         false
%
% The returned parameters structure also includes the following derived fields:
%
%   numberOfChannels
%   numberOfSites
%   injectionChannels
%
% See also WSEARCH, WEVENT, and WPROPERTIES.

% Shourov K. Chatterji <shourov@ligo.caltech.edu>

% $Id: wparameters.m 2314 2009-09-08 07:05:11Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     process command line arguments                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(0, 2, nargin));

% apply default arguments
if (nargin < 1) || isempty(parameterFile),
  parameterFile = './parameters.txt';
end
if (nargin < 2) || isempty(debugLevel),
  debugLevel = 1;
end

% validate command line arguments
if ~exist(parameterFile, 'file'),
  error('could not find parameter file:\n%s\n', parameterFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            initialize parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize parameters
channelNames = [];
frameTypes = [];

analysisMode = [];
sampleFrequency = [];
qRange = [];
frequencyRange = [];
maximumMismatch = [];
falseEventRate = [];
blockDuration = [];
conditionDuration = [];
timeShifts = [];
stateNames = [];
stateTypes = [];
stateMasks = [];
injectionNames = [];
injectionTypes = [];
injectionFactors = [];
injectionTimeShifts = [];
highPassCutoff = [];
lowPassCutoff = [];
whiteningDuration = [];
transientFactor = [];
doubleWhiten = [];
extraBlockOverlap = [];
outlierFactor = [];
maximumSignificants = [];
maximumTriggers = [];
durationInflation = [];
bandwidthInflation = [];
triggerFields = [];
triggerFormat = [];
randomSeed = [];

applyClustering = [];
clusterMethod = [];
clusterParameter1 = [];
clusterParameter2 = [];
clusterParameter3 = [];
distanceMetric = [];
writeClusters = [];

coincidenceNumber = [];
maximumCoincidents = [];

skyPosition = [];
skyCoordinateSystem = [];
applyVeto = [];
falseVetoRate = [];
uncertaintyFactor = [];
correlationFactor = [];
vetoDurationFactor = [];
vetoBandwidthFactor = [];
maximumConsistents = [];

maxFollowTriggers = [];
xCoherentCheck = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               read parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open parameter file for reading
parameterFileID = fopen(parameterFile, 'r');

% begin loop over parameter file
while ~feof(parameterFileID),

  % read one line from parameter file
  parameterLine = fgetl(parameterFileID);

  % remove any comments
  commentIndices = min([findstr(parameterLine, '#') ...
                        findstr(parameterLine, '%') ...
                        findstr(parameterLine, '//')]);
  if ~isempty(commentIndices),
    parameterLine = parameterLine(1 : (commentIndices(1) - 1));
  end

  % remove leading and trailing blanks
  parameterLine = wstrtrim(parameterLine);

  % if empty line, skip to the next line
  if isempty(parameterLine),
    continue;
  end

  % locate field separator
  colonIndex = strfind(parameterLine, ':');

  % if field separator not located, report syntax error
  if isempty(colonIndex),
    error('syntax error processing parameter file:\n%s\n', ...
          parameterLine);
  end

  % parse parameter line
  colonIndex = colonIndex(1);
  parameterName = parameterLine(1 : colonIndex - 1);
  parameterValue = parameterLine((colonIndex + 1) : end);
  parameterName = wstrtrim(parameterName);
  parameterValue = wstrtrim(parameterValue);

  % report parameter settings
  wlog(debugLevel, 1, '  %-25s%s\n', [parameterName ':'], parameterValue);

  % assign parameters based on name
  switch parameterName,

    % required parameters
    case 'channelNames',
      channelNames = eval(parameterValue);
    case 'frameTypes',
      frameTypes = eval(parameterValue);

    % optional parameters
    case 'analysisMode',
      analysisMode = eval(parameterValue);
    case 'sampleFrequency',
      sampleFrequency = eval(parameterValue);
    case 'qRange',
      qRange = eval(parameterValue);
    case 'frequencyRange',
      frequencyRange = eval(parameterValue);
    case 'maximumMismatch',
      maximumMismatch = eval(parameterValue);
    case 'falseEventRate',
      falseEventRate = eval(parameterValue);
    case 'blockDuration',
      blockDuration = eval(parameterValue);
    case 'conditionDuration',
      conditionDuration = eval(parameterValue);
    case 'timeShifts',
      timeShifts = eval(parameterValue);
    case 'stateNames',
      stateNames = eval(parameterValue);
    case 'stateTypes',
      stateTypes = eval(parameterValue);
    case 'stateMasks',
      stateMasks = eval(parameterValue);
    case 'injectionNames',
      injectionNames = eval(parameterValue);
    case 'injectionTypes',
      injectionTypes = eval(parameterValue);
    case 'injectionFactors',
      injectionFactors = eval(parameterValue);
    case 'injectionTimeShifts',
      injectionTimeShifts = eval(parameterValue);
    case 'highPassCutoff',
      highPassCutoff = eval(parameterValue);
    case 'lowPassCutoff',
      lowPassCutoff = eval(parameterValue);
    case 'whiteningDuration',
      whiteningDuration = eval(parameterValue);
    case 'transientFactor',
      transientFactor = eval(parameterValue);
    case 'doubleWhiten',
      doubleWhiten = eval(parameterValue);
    case 'extraBlockOverlap',
      extraBlockOverlap = eval(parameterValue);
    case 'outlierFactor',
      outlierFactor = eval(parameterValue);
    case 'maximumSignificants',
      maximumSignificants = eval(parameterValue);
    case 'maximumTriggers',
      maximumTriggers = eval(parameterValue);
    case 'durationInflation',
      durationInflation = eval(parameterValue);
    case 'bandwidthInflation',
      bandwidthInflation = eval(parameterValue);
    case 'triggerFields',
      triggerFields = eval(parameterValue);
    case 'triggerFormat',
      triggerFormat = eval(parameterValue);
    case 'randomSeed',
      randomSeed = eval(parameterValue);

    % clustering parameters
    case 'applyClustering',
      applyClustering = eval(parameterValue);
    case 'clusterMethod',
      clusterMethod = eval(parameterValue);
    case 'clusterRadius',
      clusterParameter1 = eval(parameterValue);
    case 'clusterDensity',
      clusterParameter2 = eval(parameterValue);
    case 'clusterSingles',
      clusterParameter3 = eval(parameterValue);
    case 'clusterLinkage',
      clusterParameter1 = eval(parameterValue);
    case 'clusterCriterion',
      clusterParameter2 = eval(parameterValue);
    case 'clusterThreshold',
      clusterParameter3 = eval(parameterValue);
    case 'distanceMetric',
      distanceMetric = eval(parameterValue);
    case 'writeClusters',
      writeClusters = eval(parameterValue);

    % 'independent' mode parameters
    case 'coincidenceNumber'
      coincidenceNumber = eval(parameterValue);
    case 'maximumCoincidents'
      maximumCoincidents = eval(parameterValue);

    % 'coherent' mode parameters
    case 'skyPosition'
      skyPosition = eval(parameterValue);
    case 'skyCoordinateSystem'
      skyCoordinateSystem = eval(parameterValue);
    case 'applyVeto',
      applyVeto = eval(parameterValue);
    case 'falseVetoRate',
      falseVetoRate = eval(parameterValue);
    case 'uncertaintyFactor',
      uncertaintyFactor = eval(parameterValue);
    case 'correlationFactor',
      correlationFactor = eval(parameterValue);
    case 'vetoDurationFactor',
      vetoDurationFactor = eval(parameterValue);
    case 'vetoBandwidthFactor',
      vetoBandwidthFactor = eval(parameterValue);
    case 'maximumConsistents',
      maximumConsistents = eval(parameterValue);

    % 'bayesian' mode parameters
    case 'maxFollowTriggers',
      maxFollowTriggers = eval(parameterValue);
    case 'xCoherentCheck',
      xCoherentCheck = eval(parameterValue);

    % handle unknown parameters
    otherwise,
     error('unknown parameter %s\n', parameterName);

  % end assign parameters based on name
  end

% end loop over parameter file entries
end

% close parameter file
fclose(parameterFileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        check for required parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test for unspecified parameters
if isempty(channelNames),
  error('channelNames not specified');
end
if isempty(frameTypes),
  error('frameTypes not specified');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            construct cell arrays                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force cell arrays and vectors
channelNames = wmat2cell(channelNames);
frameTypes = wmat2cell(frameTypes);
stateNames = wmat2cell(stateNames, ~isempty(stateNames));
stateTypes = wmat2cell(stateTypes, ~isempty(stateTypes));
injectionNames = wmat2cell(injectionNames, ~isempty(injectionNames));
injectionTypes = wmat2cell(injectionTypes, ~isempty(injectionTypes));
triggerFields = wmat2cell(triggerFields, ~isempty(triggerFields));
if iscell(timeShifts) && ~isempty(timeShifts),
  timeShifts = [timeShifts{:}];
end
if iscell(injectionFactors) && ~isempty(injectionFactors),
  injectionFactors = [injectionFactors{:}];
end
if iscell(injectionTimeShifts) && ~isempty(injectionTimeShifts),
  injectionTimeShifts = [injectionTimeShifts{:}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              set unset defaults                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of channels
numberOfChannels = length(channelNames);

% number of sites
sites = unique(regexprep(channelNames, '.:.*$', ''));
numberOfSites = length(sites);

% default values
if isempty(analysisMode),
  analysisMode='independent';
end
if isempty(sampleFrequency),
  sampleFrequency = 4096;
end
if isempty(qRange),
  qRange = [sqrt(11) 100];
end
if isempty(frequencyRange);
  frequencyRange = [48 Inf];
end
if isempty(maximumMismatch);
  maximumMismatch = 0.2;
end
if isempty(falseEventRate);
  falseEventRate = 1;
end
if isempty(blockDuration),
  blockDuration = 64;
end
if isempty(conditionDuration),
  conditionDuration = blockDuration;
end
if isempty(timeShifts),
  timeShifts = zeros(1, numberOfChannels);
end
if isempty(injectionNames),
  injectionNames = cell(1, numberOfChannels);
  [injectionNames{:}] = deal('NONE');
end
if isempty(injectionTypes),
  injectionTypes = cell(1, numberOfChannels);
  [injectionTypes{:}] = deal('NONE');
end
if isempty(injectionFactors),
  injectionFactors = zeros(1, numberOfChannels);
end
if isempty(injectionTimeShifts),
  injectionTimeShifts = zeros(1, numberOfChannels);
end
if isempty(transientFactor)
  transientFactor = 4;
end
if isempty(doubleWhiten)
  doubleWhiten = 1;
end
if isempty(extraBlockOverlap),
  extraBlockOverlap = 0;
end
if isempty(outlierFactor),
  outlierFactor = 2.0;
end
if isempty(maximumSignificants),
  maximumSignificants = 1e5;
end
if isempty(maximumTriggers),
  maximumTriggers = 1e3;
end
if isempty(durationInflation),
  durationInflation = 1.0;
end
if isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
end
if isempty(triggerFormat),
  triggerFormat = 'txt';
end
if isempty(randomSeed),
  randomSeed = sum(1e6 * clock);
end

% defaults for clustering
if isempty(applyClustering),
  applyClustering = 0;
end
if isempty(clusterMethod),
  clusterMethod = 'density';
end
switch lower(clusterMethod),
  case 'density',
    if isempty(clusterParameter1),
      clusterParameter1 = 4.0;
    end
    if isempty(clusterParameter2),
      clusterParameter2 = 3.0;
    end
    if isempty(clusterParameter3),
      clusterParameter3 = 1;
    end
  case 'hierarchical',
    if isempty(clusterParameter1),
      clusterParameter1 = 'single';
    end
    if isempty(clusterParameter2),
      clusterParameter2 = 'distance';
    end
    if isempty(clusterParameter3),
      clusterParameter3 = 4.0;
    end
end
if isempty(distanceMetric),
  distanceMetric = 'integratedMismatch';
end
if isempty(writeClusters),
  writeClusters = 0;
end

% defaults for independent analysis
if isempty(coincidenceNumber),
  coincidenceNumber = 0;
end
if isempty(maximumCoincidents),
  maximumCoincidents = Inf;
end

% defaults for coherent analysis
if isempty(skyPosition),
  skyPosition = [];
end
if isempty(skyCoordinateSystem)
  skyCoordinateSystem = 'equatorial';
end
if isempty(applyVeto),
  applyVeto = 1;
end
if isempty(falseVetoRate),
  falseVetoRate = 0;
end
if isempty(uncertaintyFactor),
  uncertaintyFactor = 0;
end
if isempty(correlationFactor),
  correlationFactor = 0;
end
if isempty(vetoDurationFactor),
  vetoDurationFactor = 0.5;
end
if isempty(vetoBandwidthFactor),
  vetoBandwidthFactor = 0.5;
end
if isempty(maximumConsistents),
  maximumConsistents = 1e3;
end
if isempty(maxFollowTriggers),
  maxFollowTriggers = 5;
end
if isempty(xCoherentCheck),
  xCoherentCheck = false;
end

% default trigger fields
if isempty(triggerFields),
  triggerFields = {'time', 'frequency', 'duration', 'bandwidth', 'normalizedEnergy'};
  if any(strcmpi(analysisMode, {'coherent'})),
    triggerFields{end + 1} = 'incoherentEnergy';
  end
  if applyClustering && ~writeClusters,
    triggerFields{end + 1} = 'clusterSize';
    triggerFields{end + 1} = 'clusterNormalizedEnergy';
  end
end

% channels with requested injections
injectionChannels = find(~strcmp(upper(injectionNames), 'NONE') & ...
                         ~strcmp(upper(injectionTypes), 'NONE') & ...
                         (injectionFactors ~= 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              reshape parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force row arrays and vectors
channelNames = channelNames(:);
frameTypes = frameTypes(:);
timeShifts = timeShifts(:);
stateNames = stateNames(:);
stateTypes = stateTypes(:);
stateMasks = stateMasks(:);
injectionNames = injectionNames(:);
injectionTypes = injectionTypes(:);
injectionFactors = injectionFactors(:);
injectionTimeShifts = injectionTimeShifts(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             validate parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate analysis mode
sites = unique(regexprep(channelNames, '.:.*$', ''));
numberOfSites = length(sites);
switch lower(analysisMode),
  case 'independent',
  case 'coherent',
    if (numberOfChannels < 2),
      error('coherent analysis mode requires at least two channels');
    end
  case 'bayesian',
    if (numberOfChannels < 2),
      error('bayesian analysis mode requires at least two channels');
    end
    if (coincidenceNumber < 2),
      error('bayesian analysis mode requires coincidenceNumber > 1');
    end
  otherwise,
    error('unknown analysis mode "%s"\n', analysisMode);
end

% validate number of frame types
if length(frameTypes) ~= numberOfChannels,
  error('number of frameTypes and channelNames are inconsistent');
end

% validate number of time shifts
if (length(timeShifts) ~= numberOfChannels) && ~isempty(timeShifts),
  error('number of timeShifts and channelNames are inconsistent');
end

% validate number of state names
if (length(stateNames) ~= numberOfChannels) && ~isempty(stateNames),
  error('number of stateNames and channelNames are inconsistent');
end

% validate number of state types
if (length(stateTypes) ~= numberOfChannels) && ~isempty(stateTypes),
  error('number of stateTypes and channelNames are inconsistent');
end

% validate number of state masks
if (length(stateMasks) ~= numberOfChannels) && ~isempty(stateMasks),
  error('number of stateMasks and channelNames are inconsistent');
end

% validate number of injection names
if (length(injectionNames) ~= numberOfChannels) && ~isempty(injectionNames),
  error('number of injectionNames and channelNames are inconsistent');
end

% validate number of injection types
if (length(injectionTypes) ~= numberOfChannels) && ~isempty(injectionTypes),
  error('number of injectionTypes and channelNames are inconsistent');
end

% validate number of injection factors
if (length(injectionFactors) ~= numberOfChannels) && ~isempty(injectionFactors),
  error('number of injectionFactors and channelNames are inconsistent');
end

% validate number of injection time shifts
if (length(injectionTimeShifts) ~= numberOfChannels) && ...
      ~isempty(injectionTimeShifts),
  error('number of injectionTimeShifts and channelNames are inconsistent');
end

% validate block duration
if rem(blockDuration,1) ~= 0,
  error('block duration must be an integer');
end

% validate sky coordinate system and position
if ~any(strcmpi(skyCoordinateSystem, {'equatorial', 'geocentric', 'galactic'})),
  error('unknown sky coordinate system');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        construct parameters structure                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% insert parameters into parameters structure
parameters.channelNames = channelNames;
parameters.frameTypes = frameTypes;

parameters.analysisMode = analysisMode;
parameters.sampleFrequency = sampleFrequency;
parameters.qRange = qRange;
parameters.frequencyRange = frequencyRange;
parameters.maximumMismatch = maximumMismatch;
parameters.falseEventRate = falseEventRate;
parameters.blockDuration = blockDuration;
parameters.conditionDuration = conditionDuration;
parameters.timeShifts = timeShifts;
parameters.stateNames = stateNames;
parameters.stateTypes = stateTypes;
parameters.stateMasks = stateMasks;
parameters.injectionNames = injectionNames;
parameters.injectionTypes = injectionTypes;
parameters.injectionFactors = injectionFactors;
parameters.injectionTimeShifts = injectionTimeShifts;
parameters.highPassCutoff = highPassCutoff;
parameters.lowPassCutoff = lowPassCutoff;
parameters.whiteningDuration = whiteningDuration;
parameters.transientFactor = transientFactor;
parameters.doubleWhiten = doubleWhiten;
parameters.extraBlockOverlap = extraBlockOverlap;
parameters.outlierFactor = outlierFactor;
parameters.falseVetoRate = falseVetoRate;
parameters.uncertaintyFactor = uncertaintyFactor;
parameters.correlationFactor = correlationFactor;
parameters.maximumSignificants = maximumSignificants;
parameters.maximumTriggers = maximumTriggers;
parameters.durationInflation = durationInflation;
parameters.bandwidthInflation = bandwidthInflation;
parameters.triggerFields = triggerFields;
parameters.triggerFormat = triggerFormat;
parameters.randomSeed = randomSeed;
parameters.numberOfChannels = numberOfChannels;
parameters.numberOfSites = numberOfSites;
parameters.injectionChannels = injectionChannels;

parameters.applyClustering = applyClustering;
parameters.clusterMethod = clusterMethod;
parameters.clusterParameter1 = clusterParameter1;
parameters.clusterParameter2 = clusterParameter2;
parameters.clusterParameter3 = clusterParameter3;
parameters.distanceMetric = distanceMetric;
parameters.writeClusters = writeClusters;

parameters.coincidenceNumber = coincidenceNumber;
parameters.maximumCoincidents = maximumCoincidents;

parameters.skyPosition = skyPosition;
parameters.skyCoordinateSystem = skyCoordinateSystem;
parameters.applyVeto = applyVeto;
parameters.falseVetoRate = falseVetoRate;
parameters.uncertaintyFactor = uncertaintyFactor;
parameters.correlationFactor = correlationFactor;
parameters.vetoDurationFactor = vetoDurationFactor;
parameters.vetoBandwidthFactor = vetoBandwidthFactor;
parameters.maximumConsistents = maximumConsistents;
parameters.maxFollowTriggers = maxFollowTriggers;
parameters.xCoherentCheck = xCoherentCheck;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return to calling function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
