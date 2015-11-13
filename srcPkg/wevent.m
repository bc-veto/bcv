function wevent(eventTime, duration, parameterFile, frameCacheFile, ...
                outputDirectory, refEventFile, debugLevel)
% WEVENT Top level function for Omega Pipeline event analysis
%
% WEVENT runs a full analysis, including followup, around the specified event
% time, and produces figures of event properties and triggers, as well as an
% html report.
%
% usage: wevent(eventTime, duration, parameterFile, frameCacheFile, ...
%               outputDirectory, debugLevel);
%
%   eventTime             gps event time
%   duration              plot duration
%   parameterFile         parameter file
%   frameCacheFile        readframedata formatted frame cache file
%   outputDirectory       directory to write results
%   refEventFile          event file in which to find event block info
%   debugLevel            verboseness level of debug output
%
% If no output directory is specified, WEVENT places the resulting trigger,
% event, image and html report files in a subdirectory of the current directory,
% named after the event center time:
%
%   events/
%     <eventTime>/
%       <channelNameX>.txt
%       <channelNameX>.png
%       ...
%       events.txt
%       skymap.txt
%       skymap.png
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
% See also WSEARCH, WTILE, WREADDATA, WRESAMPLE, WCONDITION, WTRANSFORM,
% WTHRESHOLD, WSELECT, and WWRITETRIGGERS.

% Authors:
% Shourov K. Chatterji <shourov@ligo.caltech.edu>
% Leo C. Stein <lstein@ligo.mit.edu>
% Jameson Rollins <jrollins@phys.columbia.edu>
% Antony Searle <antony.searle@anu.edu.au>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             start analysis timer                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start time of analysis
analysisTimer = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         parse command line arguments                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 7, nargin));

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
if (nargin < 6),
  refEventFile = [];
end
if (nargin < 7) || isempty(debugLevel),
  debugLevel = 1;
end

% ensure numeric event time, duration, and debug level
if ischar(eventTime),
  eventTime = str2double(eventTime);
end
if ischar(duration),
  duration = str2double(duration);
end
if ischar(debugLevel),
  debugLevel = str2double(debugLevel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 write header                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report analysis name and time
wlog(debugLevel, 1, 'Omega Event analysis\n');
wlog(debugLevel, 1, 'Run by %s on %s at %s\n', ...
     getenv('USER'), datestr(analysisTimer, 29), datestr(analysisTimer, 13));

% report command line arguments
wlog(debugLevel, 1, '  event time:              %f\n', eventTime);
wlog(debugLevel, 1, '  plot duration:           %f\n', duration);
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
durationInflation = parameters.durationInflation;
bandwidthInflation = parameters.bandwidthInflation;
triggerFormat = parameters.triggerFormat;
randomSeed = parameters.randomSeed;

% always generate a report
generateReport = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       initialize random number generators                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set random number generator seeds
rand('state', randomSeed);
randn('state', randomSeed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          determine times and report                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that the specified duration is not longer than the block duration
if duration > blockDuration,
  error('specified duration is longer than blockDuration')
end

% event string
eventString = sprintf('%#020.9f', eventTime);

% f an event file has been specified, determine the block start and stop time by
% finding the appropriate matching event in the file.  this is so that wevent
% can analyze the exact same block
if ~isempty(refEventFile),

  wlog(debugLevel, 1, 'determining block from event file:\n');
  wlog(debugLevel, 1, '  %s\n', refEventFile);

  % read in event file
  [ref.dummy, ref.blockStartTime, ref.blockStopTime, ref.time] = ...
      textread(refEventFile, '%s %f %f %f %*[^\n]');

  % find the time that most closely matches the one specified to wevent
  [dummy, index] = min(abs(eventTime - ref.time));

  % pull out the block start/stop time
  blockStartTime = ref.blockStartTime(index);
  blockStopTime = ref.blockStopTime(index);

  clear ref dummy index;

else

  % center block around event time (rounded to nearest integer)
  blockStartTime = round(eventTime - (blockDuration / 2));

  % force start time alignment with data samples
  blockStartTime = floor(blockStartTime * sampleFrequency) / ...
                   sampleFrequency;

  % stop time of block
  blockStopTime = blockStartTime + blockDuration;

end

% report status
wlog(debugLevel, 1, '  block start time:    %9.2f\n', blockStartTime);
wlog(debugLevel, 1, '  block stop time:     %9.2f\n', blockStopTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              output directory                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if outputDirectory not specified, make one based on center time
if isempty(outputDirectory),
  outputDirectory = sprintf('events/%s', eventString);
end
    
% report status
wlog(debugLevel, 1, 'creating output directory\n');
wlog(debugLevel, 1, '  outputDirectory:         %s\n', outputDirectory);

% create event directory
unix(['mkdir -p ' outputDirectory]);

% copy in parameter file
unix(['cp -f ' parameterFile ' ' outputDirectory '/parameters.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              output file paths                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time string
fileTimeString = sprintf('%#010d-%d', floor(blockStartTime), blockDuration);

% generate output file path
outputFiles = wgenfilepaths(channelNames, fileTimeString, ...
                            outputDirectory, triggerFormat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            read frame file cache                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'reading frame cache\n');

% read frame file cache
frameCache = loadframecache(frameCacheFile);

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
%                               block analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wlog(debugLevel, 1, 'block analysis\n');

[triggersThresholded, triggersDownselected, triggersVetoed, ...
 triggersClustered, triggersCoincident, event, skymap, transforms] = ...
    wblock(blockStartTime, tiling, eventTime, ...
           parameters, frameCache, outputFiles, debugLevel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          generate report with plots                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate report only if requested to do so
if generateReport,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               open html report                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open html report
htmlFile = [outputDirectory '/index.html'];
htmlFID = fopen(htmlFile, 'wt');

% begin html report
fprintf(htmlFID, '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">\n');
fprintf(htmlFID, '<html>\n');
fprintf(htmlFID, '<head>\n');
fprintf(htmlFID, '<title>Event %s</title>\n', eventString);
fprintf(htmlFID, '<link rel="icon" type="image/x-icon" href="omega.logo.icon.png" />\n');
fprintf(htmlFID, '<link rel="stylesheet" type="text/css" href="wstyle.css" />\n');
fprintf(htmlFID, '</head>\n');
fprintf(htmlFID, '<body>\n');
fprintf(htmlFID, '\n');

fprintf(htmlFID, '<h1 class="title"><a href="https://geco.phys.columbia.edu/omega" class="title">Omega Pipeline</a> Event</h1>\n');
fprintf(htmlFID, '<h1>Event %s</h1>\n', eventString);
fprintf(htmlFID, '\n');
fprintf(htmlFID, '<div class="main">\n');
fprintf(htmlFID, '\n');

fprintf(htmlFID, '<div class="section">\n');
fprintf(htmlFID, '<ul>\n');
fprintf(htmlFID, '<li><a href="parameters.txt">parameters.txt</a></li>\n');
fprintf(htmlFID, '<li><a href="log.txt">log.txt</a></li>\n');
fprintf(htmlFID, '</ul>\n');
fprintf(htmlFID, '</div>\n');
fprintf(htmlFID, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          setup figure properties                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print all figures to files
printFigs = true;
% use version 0 printing
printVers = 0;
figureHandle = figure;
%figureHandle = figure('visible', 'off');

% {spectro,event}gram parameters
referenceTime = eventTime;
timeRange = [-1 +1] * duration / 2;
normalizedEnergyRange = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   produce bayesian skymap and event info                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(skymap)

  % initialize html block
  fprintf(htmlFID, '<h2>Bayesian followup</h2>\n');
  fprintf(htmlFID, '<div class="section">\n');
  fprintf(htmlFID, '<div class="subsection">\n');
  fprintf(htmlFID, '<table><tr>\n');

  % report status
  wlog(debugLevel, 1, '  plotting bayesian skymap\n');

  % plot skymap
  close all;
  figureHandle = figure;
  whybridskyplot(skymap);

  % print skymap
  [figPath, thumbPath] = ...
      wprintfig(figureHandle, ...
                [outputDirectory '/skymap'], ...
                printFigs, printVers);

  % write skymap html
  fprintf(htmlFID, '<td>\n');
  fprintf(htmlFID, '<h3><a href="%s">Sky map</a></h3>\n', ...
          strrep(basename(outputFiles.SKYMAP_BASE), ...
                 '@TIME@', sprintf('%#020.9f', event.time)));

  fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a>\n', ...
          basename(figPath), basename(thumbPath));
  fprintf(htmlFID, '</td>\n');

  % write event properties html
  fprintf(htmlFID, '<td>\n');
  fprintf(htmlFID, '<h3><a href="%s">Event info</a></h3>\n', ...
          basename(outputFiles.EVENTS));
  fprintf(htmlFID, '<table>\n');
  fprintf(htmlFID, '<tr><td><b>probSignal:</b></td><td><b>%f</b></td></tr>\n', event.probSignal);
  fprintf(htmlFID, '<tr><td>probGlitch:</td><td>%f</td></tr>\n', event.probGlitch);
  fprintf(htmlFID, '<tr><td colspan="2" align="center"><hr /></td></tr>\n');
  fprintf(htmlFID, '<tr><td>frequency:</td><td>%f</td></tr>\n', event.frequency);
  fprintf(htmlFID, '<tr><td>duration:</td><td>%f</td></tr>\n', event.duration);
  fprintf(htmlFID, '<tr><td>bandwidth:</td><td>%f</td></tr>\n', event.bandwidth);
  fprintf(htmlFID, '<tr><td colspan="2" align="center"><hr /></td></tr>\n');
  fprintf(htmlFID, '<tr><td>modeTheta:</td><td>%f</td></tr>\n', event.modeTheta);
  fprintf(htmlFID, '<tr><td>modePhi:</td><td>%f</td></tr>\n', event.modePhi);

  fprintf(htmlFID, '</table>\n');
  fprintf(htmlFID, '</td>\n');

  % finalize html block
  fprintf(htmlFID, '</tr></table>\n');
  fprintf(htmlFID, '</div>\n'); % end section
  fprintf(htmlFID, '</div>\n'); % end subsection
  fprintf(htmlFID, '\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             produce spectrograms                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  plotting spectrograms\n');

% intialize html block
fprintf(htmlFID, '<h2>Spectrograms of transformed data</h2>\n');
fprintf(htmlFID, '<div class="section">\n');
fprintf(htmlFID, '<div class="subsection">\n');
fprintf(htmlFID, '<table><tr>\n');

% begin loop over channels
for channelNumber = 1 : length(transforms),

  % initialize html block for this channel
  fprintf(htmlFID, '<td>\n');
  fprintf(htmlFID, '<h3>%s</h3>\n', ...
          transforms{channelNumber}.channelName);
  figNameBase = [transforms{channelNumber}.channelName '_spectrogram_'];

  % begin loop over q planes
  for plane = 1 : tiling.numberOfPlanes,

    % plot spectrogram
    wspectrogram(transforms{channelNumber}, tiling, blockStartTime, ...
                 referenceTime, timeRange, frequencyRange, ...
                 tiling.planes{plane}.q, normalizedEnergyRange);

    % print spectrogram
    [figPath, thumbPath] = ...
        wprintfig(figureHandle, ...
                  [outputDirectory '/' figNameBase int2str(plane)], ...
                  printFigs,printVers);

    % add spectrogram to html report
    fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a><br />\n', ...
            basename(figPath), basename(thumbPath));

  % end loop over q planes
  end

  fprintf(htmlFID, '</td>\n');
  
% end loop over channels
end

% finalize html block
fprintf(htmlFID, '</tr></table>\n');
fprintf(htmlFID, '</div>\n'); % end section
fprintf(htmlFID, '</div>\n'); % end subsection
fprintf(htmlFID, '\n');

% clear the transforms as they are only needed for the spectrogram
clear transforms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        produce thresholded eventgrams                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  plotting thresholded eventgrams\n');

% initialize html block
fprintf(htmlFID, '<h2>Eventgrams of significant tiles</h2>\n');
fprintf(htmlFID, '<div class="section">\n');
fprintf(htmlFID, '<div class="subsection">\n');
fprintf(htmlFID, '<table><tr>\n');

% begin loop over channels
for channelNumber = 1 : length(triggersThresholded),

  % plot eventgram
  weventgram(triggersThresholded{channelNumber}, tiling, blockStartTime, ...
             referenceTime, ...
             timeRange, frequencyRange, durationInflation, ...
             bandwidthInflation, normalizedEnergyRange);

  % print eventgram
  [figPath, thumbPath] = ...
      wprintfig(figureHandle, ...
                [outputDirectory '/' ...
                 triggersThresholded{channelNumber}.channelName], ...
                printFigs,printVers);

  % write eventgram html
  fprintf(htmlFID, '<td>\n');
  fprintf(htmlFID, '<h3>%s</h3>\n', ...
          triggersThresholded{channelNumber}.channelName);
  fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a><br />\n', ...
          basename(figPath), basename(thumbPath));
  fprintf(htmlFID, '</td>\n');

% end loop over channels
end

% finalize html block
fprintf(htmlFID, '</tr></table>\n');
fprintf(htmlFID, '</div>\n'); % end section
fprintf(htmlFID, '</div>\n'); % end subsection
fprintf(htmlFID, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      produce downselected eventgrams                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, '  plotting downselected eventgrams\n');

% initialize html block
fprintf(htmlFID, '<h2>Eventgrams of non-overlapping significant tiles</h2>\n');
fprintf(htmlFID, '<div class="section">\n');
fprintf(htmlFID, '<div class="subsection">\n');
fprintf(htmlFID, '<table><tr>\n');

% begin loop over channels
for channelNumber = 1 : length(triggersDownselected),

  % plot eventgram
  weventgram(triggersDownselected{channelNumber}, tiling, blockStartTime, ...
             referenceTime, ...
             timeRange, frequencyRange, durationInflation, ...
             bandwidthInflation, normalizedEnergyRange);

  % print eventgram
  [figPath, thumbPath] = ...
      wprintfig(figureHandle, ...
                [outputDirectory '/' ...
                 triggersDownselected{channelNumber}.channelName], ...
                printFigs,printVers);

  % write eventrgram html
  fprintf(htmlFID, '<td>\n');
  fprintf(htmlFID, '<h3><a href="%s">%s</a></h3>\n', ...
          basename(outputFiles.DOWNSELECT{channelNumber}), ...
          triggersDownselected{channelNumber}.channelName);
  fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a>\n', ...
          basename(figPath), basename(thumbPath));
  fprintf(htmlFID, '</td>\n');

% end loop over channels
end

% finalize html block
fprintf(htmlFID, '</tr></table>\n');
fprintf(htmlFID, '</div>\n'); % end section
fprintf(htmlFID, '</div>\n'); % end subsection
fprintf(htmlFID, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           produce vetoed eventgrams                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(triggersVetoed),

  % report status
  wlog(debugLevel, 1, '  plotting veto eventgrams\n');

  % initialize html block
  fprintf(htmlFID, '<h2>Eventgrams of non-overlapping consistent tiles</h2>\n');
  fprintf(htmlFID, '<div class="section">\n');
  fprintf(htmlFID, '<div class="subsection">\n');
  fprintf(htmlFID, '<table><tr>\n');

  % begin loop over channels
  for channelNumber = 1 : length(triggersVetoed),

    % plot eventgram
    weventgram(triggersVetoed{channelNumber}, tiling, blockStartTime, ...
               referenceTime, ...
               timeRange, frequencyRange, durationInflation, ...
               bandwidthInflation, normalizedEnergyRange);

    % print eventgram
    [figPath, thumbPath] = ...
        wprintfig(figureHandle, ...
                  [outputDirectory '/' ...
                   triggersVetoed{channelNumber}.channelName], ...
                  printFigs,printVers);

    % write eventgram html
    fprintf(htmlFID, '<td>\n');
    fprintf(htmlFID, '<h3>%s</h3>\n', ...
            triggersVetoed{channelNumber}.channelName);
    fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a>\n', ...
            basename(figPath), basename(thumbPath));
    fprintf(htmlFID, '</td>\n');

  % end loop over channels
  end

  % finalize html block
  fprintf(htmlFID, '</tr></table>\n');
  fprintf(htmlFID, '</div>\n'); % end section
  fprintf(htmlFID, '</div>\n'); % end subsection
  fprintf(htmlFID, '\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          produce cluster eventgrams                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(triggersClustered),

  % report status
  wlog(debugLevel, 1, '  plotting cluster eventgrams\n');

  % initialize html block
  fprintf(htmlFID, '<h2>Eventgrams of clustered tiles</h2>\n');
  fprintf(htmlFID, '<div class="section">\n');
  fprintf(htmlFID, '<div class="subsection">\n');
  fprintf(htmlFID, '<table><tr>\n');

  % begin loop over channels
  for channelNumber = 1 : length(triggersClustered),

    % plot eventgram
    weventgram(triggersClustered{channelNumber}, tiling, blockStartTime, ...
               referenceTime, ...
               timeRange, frequencyRange, durationInflation, ...
               bandwidthInflation, normalizedEnergyRange);

    % print eventgram
    [figPath, thumbPath] = ...
        wprintfig(figureHandle, ...
                  [outputDirectory '/' ...
                   triggersClustered{channelNumber}.channelName], ...
                  printFigs,printVers);

    % write eventgram html
    fprintf(htmlFID, '<td>\n');
    fprintf(htmlFID, '<h3><a href="%s">%s</a></h3>\n', ...
            basename(outputFiles.CLUSTER{channelNumber}), ...
            triggersClustered{channelNumber}.channelName);
    fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a>\n', ...
            basename(figPath), basename(thumbPath));
    fprintf(htmlFID, '</td>\n');

  % end loop over channels
  end

  % finalize html block
  fprintf(htmlFID, '</tr></table>\n');
  fprintf(htmlFID, '</div>\n'); % end section
  fprintf(htmlFID, '</div>\n'); % end subsection
  fprintf(htmlFID, '\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       produce coincident eventgrams                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(triggersCoincident),

  % report status
  wlog(debugLevel, 1, '  plotting coincident eventgrams\n');

  % initialize html block
  fprintf(htmlFID, '<h2>Eventgrams of coincident tiles</h2>\n');
  fprintf(htmlFID, '<div class="section">\n');
  fprintf(htmlFID, '<div class="subsection">\n');
  fprintf(htmlFID, '<table><tr>\n');
  
  % begin loop over channels
  for channelNumber = 1 : length(triggersCoincident),

    % plot eventgram
    weventgram(triggersCoincident{channelNumber}, tiling, blockStartTime, ...
               referenceTime, ...
               timeRange, frequencyRange, durationInflation, ...
               bandwidthInflation, normalizedEnergyRange);

    % print eventgram
    [figPath, thumbPath] = ...
        wprintfig(figureHandle, ...
                  [outputDirectory '/' ...
                   triggersCoincident{channelNumber}.channelName], ...
                  printFigs,printVers);

    % write eventgram html
    fprintf(htmlFID, '<td>\n');
    fprintf(htmlFID, '<h3><a href="%s">%s</a></h3>\n', ...
            basename(outputFiles.COINCIDE{channelNumber}), ...
            triggersCoincident{channelNumber}.channelName);
    fprintf(htmlFID, '<a href="./%s"><img src="./%s"></a>\n', ...
            basename(figPath), basename(thumbPath));
    fprintf(htmlFID, '</td>\n');

  % end loop over channels
  end

  % finalize html block
  fprintf(htmlFID, '</tr></table>\n');
  fprintf(htmlFID, '</div>\n'); % end section
  fprintf(htmlFID, '</div>\n'); % end subsection
  fprintf(htmlFID, '\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              close html report                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end html report
fprintf(htmlFID, '</div>\n');
fprintf(htmlFID, '\n');
fprintf(htmlFID, '<hr>\n');
fprintf(htmlFID, '\n');
fprintf(htmlFID, '<div class="footer">\n');
fprintf(htmlFID, 'Created by user %s on %s at %s<br />\n', ...
        getenv('USER'), datestr(clock, 29), datestr(clock, 13));
fprintf(htmlFID, '</div>\n');
fprintf(htmlFID, '\n');
fprintf(htmlFID, '</body>\n');
fprintf(htmlFID, '</html>\n');

% close html report
fclose(htmlFID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              end generate report                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 write footer                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
wlog(debugLevel, 1, 'event analysis complete\n');

% report total elapsed time
wlog(debugLevel, 1, 'elapsed time:            %5.0f seconds\n', ...
        etime(clock, analysisTimer));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure all figures are closed
close all;

% return to calling function
return;
