function wsearch(startTime, stopTime, parameterFile, frameCacheFile, ...
                 outputDirectory, debugLevel)
% WSEARCH Top level function for the Omega Pipeline burst search
%
% WSEARCH applies the discrete Q transform to search for statistically
% significant transient events in data from interferometric
% gravitational-wave detectors.
%
% usage: wsearch(startTime, stopTime, parameterFile, frameCacheFile, ...
%                outputDirectory, debugLevel)
%
%   startTime           gps start time of analysis
%   stopTime            gps stop time of analysis
%   parameterFile       parameter file
%   frameCacheFile      readframedata formatted frame cache file
%   outputDirectory     directory to write results
%   debugLevel          verboseness level of debug output
%
% If the specified stop time is less than the specified start time, it is
% instead assumed to be the duration of the analysis.  Non-integer start
% and stop times are truncated to the nearest integer.
%
% If no output directory is specified, WSEARCH places the resulting trigger and
% event files in a subdirectory of the current directory, named after the
% specified segment:
%
%   segments/
%     <startTime>-<stopTime>/
%       livetime.txt
%       <channelName1>.txt
%       <channelName2>.txt
%       ...
%
% If an outputDirectory is specified, then all output is written to that
% directory instead.
%
% If no parameter file or frame cache file is specified, WSEARCH looks
% for the files parameters.txt or framecache.txt in the current directory.
%
% The specified debugLevel controls the amount of detail in the output log.
% A debugLevel of unity is assumed by default.
%
% See also WBLOCK, WEVENT

% Authors:
% Shourov K. Chatterji <shourov@ligo.caltech.edu>
% Leo C. Stein <lstein@ligo.mit.edu>
% Jameson Graef Rollins <jrollins@phys.columbia.edu>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             start analysis timer                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start time of analysis
analysisTimer = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         parse command line arguments                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 6, nargin));

% apply default arguments
if (nargin < 3) || isempty(parameterFile),
  parameterFile = 'parameters.txt';
end
if (nargin < 4) || isempty(frameCacheFile),
  frameCacheFile = 'framecache.txt';
end
if (nargin < 5),
  outputDirectory = [];
end
if (nargin < 6) || isempty(debugLevel),
  debugLevel = 1;
end

% ensure numeric start time, stop time, and debug level
if ischar(startTime),
  startTime = str2double(startTime);
end
if ischar(stopTime),
  stopTime = str2double(stopTime);
end
if ischar(debugLevel),
  debugLevel = str2double(debugLevel);
end

% truncate to integer start and stop times
startTime = floor(startTime);
stopTime = floor(stopTime);

% handle duration argument in place of stop time
if stopTime < startTime,
  stopTime = startTime + stopTime;
end

% segment duration
segmentDuration = stopTime - startTime;

% string start and stop times
startTimeString = sprintf('%#010d', startTime);
stopTimeString = sprintf('%#010d', stopTime);
segmentDurationString = sprintf('%d', segmentDuration);

% creat output path if not specified
if isempty(outputDirectory),
    % path to segment output directory
    outputDirectory = sprintf('./segments/%s-%s', startTimeString, stopTimeString);
% check if specified path is absolute
elseif ~strcmp(outputDirectory(1),'/'),
    % make explicitly relative to current directory if not absolute
    outputDirectory = ['./' outputDirectory];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 write header                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report analysis name and time
wlog(debugLevel, 1, 'Omega Search analysis\n');
wlog(debugLevel, 1, 'Run by %s on %s at %s\n', ...
     getenv('USER'), datestr(analysisTimer, 29), datestr(analysisTimer, 13));

% report command line arguments
wlog(debugLevel, 1, '  analysis start time:     %s\n', startTimeString);
wlog(debugLevel, 1, '  analysis stop time:      %s\n', stopTimeString);
wlog(debugLevel, 1, '  parameter file:          %s\n', parameterFile);
wlog(debugLevel, 1, '  framecache file:         %s\n', frameCacheFile);
wlog(debugLevel, 1, '  output directory:        %s\n', outputDirectory);
wlog(debugLevel, 1, '  debug level:             %d\n', debugLevel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               read parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'reading parameter file\n');

% read parameters
parameters = wparameters(parameterFile, debugLevel);

% extract needed parameters from parameters structure
channelNames = parameters.channelNames;
blockDuration = parameters.blockDuration;
sampleFrequency = parameters.sampleFrequency;
frequencyRange = parameters.frequencyRange;
qRange = parameters.qRange;
maximumMismatch = parameters.maximumMismatch;
highPassCutoff = parameters.highPassCutoff;
lowPassCutoff = parameters.lowPassCutoff;
whiteningDuration = parameters.whiteningDuration;
transientFactor = parameters.transientFactor;
timeShifts = parameters.timeShifts;
extraBlockOverlap = parameters.extraBlockOverlap;
triggerFormat = parameters.triggerFormat;
randomSeed = parameters.randomSeed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       initialize random number generators                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set random number generator seeds
rand('state', randomSeed);
randn('state', randomSeed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              tile signal space                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'tiling search space\n');

% generate Q transform tiling
tiling = wtile(blockDuration, qRange, frequencyRange, sampleFrequency, ...
               maximumMismatch, highPassCutoff, lowPassCutoff, ...
               whiteningDuration, transientFactor);

wlog(debugLevel, 2, '  whiteningDuration:   %9.2f\n', tiling.whiteningDuration);
wlog(debugLevel, 2, '  transientDuration:   %9.2f\n', tiling.transientDuration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              partition segment                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'portioning segment\n');

% force time shift alignment with sample interval
timeShifts = round(timeShifts * sampleFrequency) / sampleFrequency;

% maximum and minimum time shifts
maximumTimeShift = max(timeShifts);
minimumTimeShift = min(timeShifts);

% determine minimum block overlap
minimumBlockOverlap = 2 * tiling.transientDuration + extraBlockOverlap;

wlog(debugLevel, 2, '  minimumBlockOverlap: %9.2f\n', minimumBlockOverlap);

% total segment livetime loss due to time shifts
timeShiftLoss = maximumTimeShift - minimumTimeShift;

% if segment duration is less than block duration
if (segmentDuration - timeShiftLoss) < blockDuration,

  % write error to log file
  error('segment too short to analyze');

end

% number of blocks in segments
numberOfBlocks = ceil((segmentDuration - timeShiftLoss - ...
                       minimumBlockOverlap) / ...
                      (blockDuration - minimumBlockOverlap));

% if only one block,
if numberOfBlocks == 1,

  % ignore block overlap
  blockOverlap = 0;

% otherwise,
else

  % overlap in seconds between blocks
  blockOverlap = (segmentDuration - timeShiftLoss - ...
                  blockDuration * numberOfBlocks) ./ ...
                 (1 - numberOfBlocks);
end

% initial blocks to analyze
blocksToAnalyze = 1 : numberOfBlocks;

% initial analyzed blocks
analyzedBlocks = [];

% name of livetime file
livetimeFile = [outputDirectory '/livetime.txt'];

% report segment information
wlog(debugLevel, 1, '  segment duration:    %9.2f seconds\n', segmentDuration);
wlog(debugLevel, 1, '  time shift loss:     %9.2f seconds\n', timeShiftLoss);
wlog(debugLevel, 1, '  block duration:      %9.2f seconds\n', blockDuration);
wlog(debugLevel, 1, '  block overlap:       %9.2f seconds\n', blockOverlap);
wlog(debugLevel, 1, '  number of blocks:    %9u\n', numberOfBlocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           create results directory                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if output directory exists
if ~exist(outputDirectory,'dir'),

  % report status
  wlog(debugLevel, 1, 'creating output directory\n');

  % create output directory
  unix(['mkdir -p ' outputDirectory]);

else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          check for analyzed blocks                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % report status
  wlog(debugLevel, 1, 'identifying previously analyzed blocks\n');

  % read list of previously analyzed blocks
  if exist(livetimeFile,'file') == 2,
    analyzedBlocks = textread(livetimeFile, '%u %*[^\n]');
  end

  % list of blocks to analyze
  blocksToAnalyze = setdiff(blocksToAnalyze, analyzedBlocks);

  % report number of previously analyzed blocks
  wlog(debugLevel, 1, '  analyzed blocks:     %9u blocks\n', length(analyzedBlocks));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              output file paths                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time string
fileTimeString = sprintf('%s-%s', startTimeString, segmentDurationString);

% generate output file path
outputFiles = wgenfilepaths(channelNames, fileTimeString, ...
                            outputDirectory, triggerFormat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            read frame file cache                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(blocksToAnalyze),

  % report status
  wlog(debugLevel, 1, 'reading frame cache\n');

  % read frame file cache
  frameCache = loadframecache(frameCacheFile);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              begin block loop                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize error flag
errorFlag = 0;

% initialize coordinate variable
coordinate = [];

% begin loop over blocks
for blockNumber = blocksToAnalyze,

  % start time of block analysis
  blockTimer = clock;

  % start time of block
  blockStartTime = startTime + maximumTimeShift + ...
                   (blockNumber - 1) * (blockDuration - blockOverlap);

  % force start time alignment with data samples
  blockStartTime = floor(blockStartTime * sampleFrequency) / ...
                   sampleFrequency;

  % stop time of block
  blockStopTime = blockStartTime + blockDuration;

  % report status
  wlog(debugLevel, 1, 'analyzing block %u of %u\n', blockNumber, numberOfBlocks);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              block analysis                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % try block and catch any errors
  try

      channelNames = ...
          wblock(blockStartTime, tiling, [], ...
                 parameters, frameCache, outputFiles, debugLevel);
  
  catch

      exception = lasterror;
      wlog(debugLevel, 1, 'ERROR: %s\n', exception.message);
      for stackNumber = 1:length(exception.stack)
        wlog(debugLevel, 1, '  in %s:%s\n', ...
             exception.stack(stackNumber).name, ...
             num2str(exception.stack(stackNumber).line));
      end

      % note error, display, and reset
      errorFlag = errorFlag + 1;
      lasterror('reset');

      % continue to next block if any errors found
      continue;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              write livetime                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % report status
  wlog(debugLevel, 1, '  writing livetime\n');

  % open livetime file for appending
  livetimeFileFID = fopen(livetimeFile, 'at');

  % test for error
  if livetimeFileFID == -1,
    error('cannot open livetime file for writing');
  end

  % analyzed channels string
  analyzedChannels = [];
  for channelNumber = 1 : length(channelNames),
      if channelNumber == 1,
          analyzedChannels = sprintf('%s',channelNames{channelNumber});
      else
          analyzedChannels = sprintf('%s,%s',analyzedChannels,channelNames{channelNumber});
      end
  end

  % write block livetime to livetime file
  fprintf(livetimeFileFID, '%#04u %#020.9f %#020.9f %s\n', ...
          blockNumber, ...
          blockStartTime + tiling.transientDuration, ...
          blockStopTime - tiling.transientDuration, ...
          analyzedChannels);

  % close livetime file
  fclose(livetimeFileFID);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              end block loop                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % report status
  wlog(debugLevel, 1, '  block complete\n');

  % report total elapsed time
  wlog(debugLevel, 1, '    elapsed time:          %5.0f seconds\n', ...
          etime(clock, blockTimer));

% end loop over blocks
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 write footer                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'segment complete\n');

% report total elapsed time
wlog(debugLevel, 1, '  elapsed time:            %5.0f seconds\n', ...
        etime(clock, analysisTimer));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exit with error status
if errorFlag,
  error('analysis had errors.');
else
  return;
end
