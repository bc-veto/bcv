function wproperties(parameterFile, triggersFile, clustersFile, ...
                     segmentsFile, injectionsFile, livetimeFile, ...
                     startTime, stopTime, channelName, timeStamp, ...
                     outputDirectory, debugLevel)
% WPROPERTIES display and print plots of trigger properties
%
% WPROPERTIES reads the specified WSEARCH unclustered and clustered
% trigger files and produces plots of the followng trigger properties
% for both the clustered and unclustered triggers
%
%   histogram of average trigger rate per minute
%   scatter plot of trigger center frequency vs. center time
%   scatter plot of trigger snr vs. center time
%   scatter plot of trigger snr vs. center frequency
%   relative distribution of triggers in frequency and Q
%   distribution of trigger rate vs. snr threshold
%   fourier transform of trigger autocorrelogram
%
% usage:
%
%   wproperties(parameterFile, triggersFile, clustersFile, ...
%               segmentsFile, injectionsFile, livetimeFile, ...
%               startTime, stopTime, channelName, timeStamp, ...
%               outputDirectory, debugLevel);
%
%   parameterFile    parameter file
%   triggersFile     unclustered triggers file
%   clustersFile     clustered triggers file
%   segmentsFile     science mode segment file
%   injectionsFile   injection mode segment file
%   livetimeFile     analyzed segment file
%   startTime        start time for computing and displaying properties
%   stopTime         stop time for computing and displaying properties
%   channelName      channel name
%   timeStamp        time stamp corresponding to start time
%   outputDirectory  output directory
%   debugLevel       verboseness level of debug output
%
% If they are not specified, the start and stop times are set to the minimum and
% maximum times listed in the livetime segments file, and the output directory
% is set to the sub-directory "properties" of the current working directory.
%
% The livetime file must consist of non-overlapping livetime segments
% sorted by increasing start time.  The livetime file should be a two
% column ASCII text list of start and stop times.
%
% The triggers and clusters files are 5 column ASCII lists of the following
% trigger properties.
%
%   central time [gps seconds]
%   central frequency [Hz]
%   duration [seconds]
%   bandwidth [Hz]
%   normalized energy []
%
% Additional columns are permitted, and will simply be ignored.
%
% See also WSEARCH.
%
% Shourov K. Chatterji <shourov@ligo.caltech.edu>
% Jameson Graef Rollins <jrollins@phys.columbia.edu>
%
% $Id: wproperties.m 1610 2009-03-29 21:29:35Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            hard coded parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% minimum signal to noise ratio threshold
rho0 = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(6, 12, nargin));

% apply default arguments
if (nargin < 7) || isempty(startTime),
  startTime = [];
end
if (nargin < 8) || isempty(stopTime),
  stopTime = [];
end
if (nargin < 9) || isempty(channelName),
  channelName = '';
end
if (nargin < 10) || isempty(timeStamp),
  timeStamp = [];
end
if (nargin < 11) || isempty(outputDirectory),
  outputDirectory = 'properties';
end
if (nargin < 12) || isempty(debugLevel),
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               read parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('reading parameter file...');

% read parameters
parameters = wparameters(parameterFile, debugLevel);

% extract parameters from parameters structure
channelNames = parameters.channelNames;
frameTypes = parameters.frameTypes;

analysisMode = parameters.analysisMode;
sampleFrequency = parameters.sampleFrequency;
blockDuration = parameters.blockDuration;
qRange = parameters.qRange;
frequencyRange = parameters.frequencyRange;
maximumMismatch = parameters.maximumMismatch;
falseEventRate = parameters.falseEventRate;
timeShifts = parameters.timeShifts;
stateNames = parameters.stateNames;
stateTypes = parameters.stateTypes;
stateMasks = parameters.stateMasks;
injectionNames = parameters.injectionNames;
injectionTypes = parameters.injectionTypes;
injectionFactors = parameters.injectionFactors;
injectionTimeShifts = parameters.injectionTimeShifts;
highPassCutoff = parameters.highPassCutoff;
lowPassCutoff = parameters.lowPassCutoff;
whiteningDuration = parameters.whiteningDuration;
doubleWhiten = parameters.doubleWhiten;
extraBlockOverlap = parameters.extraBlockOverlap;
outlierFactor = parameters.outlierFactor;
maximumSignificants = parameters.maximumSignificants;
maximumTriggers = parameters.maximumTriggers;
durationInflation = parameters.durationInflation;
bandwidthInflation = parameters.bandwidthInflation;
triggerFields = parameters.triggerFields;
triggerFormat = parameters.triggerFormat;
randomSeed = parameters.randomSeed;
numberOfChannels = parameters.numberOfChannels;
numberOfSites = parameters.numberOfSites;
injectionChannels = parameters.injectionChannels;

applyClustering = parameters.applyClustering;
clusterMethod = parameters.clusterMethod;
clusterParameter1 = parameters.clusterParameter1;
clusterParameter2 = parameters.clusterParameter2;
clusterParameter3 = parameters.clusterParameter3;
distanceMetric = parameters.distanceMetric;
writeClusters = parameters.writeClusters;

coincidenceNumber = parameters.coincidenceNumber;
maximumCoincidents = parameters.maximumCoincidents;

skyPosition = parameters.skyPosition;
skyCoordinateSystem = parameters.skyCoordinateSystem;
applyVeto = parameters.applyVeto;
falseVetoRate = parameters.falseVetoRate;
uncertaintyFactor = parameters.uncertaintyFactor;
correlationFactor = parameters.correlationFactor;
vetoDurationFactor = parameters.vetoDurationFactor;
vetoBandwidthFactor = parameters.vetoBandwidthFactor;
maximumConsistents = parameters.maximumConsistents;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                load segments                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load science segments
disp('loading science segments...');
segments = dataread('file', segmentsFile, '', 'commentstyle', 'shell');

% load injection segments
disp('loading injection segments...');
injections = dataread('file', injectionsFile, '', 'commentstyle', 'shell');

% load livetime segments
disp('loading livetime segments...');
livetime = dataread('file', livetimeFile, '', 'commentstyle', 'shell');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                load triggers                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load unclustered triggers
disp('loading unclustered triggers...');
[triggersTime, triggersFrequency, triggersDuration, ...
 triggersBandwidth, triggersNormalizedEnergy] = ...
    dataread('file', triggersFile, '%f %f %f %f %f%*[^\n]', ...
             'commentstyle', 'matlab');
triggersQ = 2 * sqrt(pi) * triggersFrequency ./ triggersBandwidth;
triggersRho = sqrt(2 * triggersNormalizedEnergy);
clear triggersDuration;
clear triggersBandwidth;
clear triggersNormalizedEnergy;

% load clustered triggers
disp('loading clustered triggers...');
[clustersTime, clustersFrequency, clustersDuration, ...
 clustersBandwidth, clustersNormalizedEnergy] = ...
    dataread('file', clustersFile, '%f %f %f %f %f%*[^\n]', ...
             'commentstyle', 'matlab');
clustersQ = 2 * sqrt(pi) * clustersFrequency ./ clustersBandwidth;
clustersRho = sqrt(2 * clustersNormalizedEnergy);
clear clustersDuration;
clear clustersBandwidth;
clear clustersNormalizedEnergy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           extract requested times                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine default start and stop times
if isempty(startTime),
  startTime = min(segments(:, 1));
end
if isempty(stopTime),
  stopTime = max(segments(:, 2));
end

% convert string times to numbers
if isstr(startTime),
  startTimeString = startTime;
  startTime = str2num(startTimeString);
else
  startTimeString = num2str(startTime);
end
if isstr(stopTime),
  stopTimeString = stopTime;
  stopTime = str2num(stopTimeString);
else
  stopTimeString = num2str(stopTime);
end

% restrict science segments to requested time period
disp('masking science segments...');
if ~isempty(segments),
  keepIndices = find((segments(:, 1) < stopTime) & ...
                     (segments(:, 2) > startTime));
  segments(:, 1) = max(startTime, segments(:, 1));
  segments(:, 2) = min(stopTime, segments(:, 2));
  segments = segments(keepIndices, :);
end

% restrict injection segments to requested time period
disp('masking injection segments...');
if ~isempty(injections),
  keepIndices = find((injections(:, 1) < stopTime) & ...
                     (injections(:, 2) > startTime));
  injections(:, 1) = max(startTime, injections(:, 1));
  injections(:, 2) = min(stopTime, injections(:, 2));
  injections = injections(keepIndices, :);
end

% restrict analyzed segments to requested time period
disp('masking livetime segments...');
if ~isempty(livetime),
  keepIndices = find((livetime(:, 1) < stopTime) & ...
                     (livetime(:, 2) > startTime));
  livetime(:, 1) = max(startTime, livetime(:, 1));
  livetime(:, 2) = min(stopTime, livetime(:, 2));
  livetime = livetime(keepIndices, :);
end

% construct non science segments
disp('constructing downtime segments');
downtimes = segments.';
downtimes = [startTime; downtimes(:); stopTime;];
downtimes = reshape(downtimes, 2, size(segments, 1) + 1).';
keepIndices = find(diff(downtimes, 1, 2) ~= 0);
downtimes = downtimes(keepIndices, :);

% number of science segments
numberOfSegments = size(segments, 1);

% number of injections segments
numberOfInjections = size(injections, 1);

% number of analyzed segments
numberOfLivetimes = size(livetime, 1);

% number of downtime segments
numberOfDowntimes = size(downtimes, 1);

% total science mode time
scienceTime = max(1, sum(diff(segments, 1, 2)));

% total injection mode time
injectionTime = max(1, sum(diff(injections, 1, 2)));

% total analyzed livetime
analysisLivetime = max(1, sum(diff(livetime, 1, 2)));

% restrict unclustered triggers to requested time period
disp('masking unclustered triggers...');
keepIndices = find((triggersTime >= startTime) & ...
                   (triggersTime < stopTime));
triggersTime = triggersTime(keepIndices);
triggersFrequency = triggersFrequency(keepIndices);
triggersQ = triggersQ(keepIndices);
triggersRho = triggersRho(keepIndices);

% restrict clustered triggers to requested time period
disp('masking clustered triggers...');
keepIndices = find((clustersTime >= startTime) & ...
                   (clustersTime < stopTime));
clustersTime = clustersTime(keepIndices);
clustersFrequency = clustersFrequency(keepIndices);
clustersQ = clustersQ(keepIndices);
clustersRho = clustersRho(keepIndices);

% renormalize time scale
disp('renormalizing time scale...');
analysisDuration = stopTime - startTime;
if analysisDuration >= 3 * 24 * 60 * 60,
  timeScale = 24 * 60 * 60;
  timeLabel = 'days';
elseif analysisDuration >= 3 * 60 * 60,
  timeScale = 60 * 60;
  timeLabel = 'hours';
elseif analysisDuration >= 3 * 60,
  timeScale = 60;
  timeLabel = 'minutes';
else
  timeScale = 1;
  timeLabel = 'seconds';
end
triggersTime = (triggersTime - startTime) / timeScale;
clustersTime = (clustersTime - startTime) / timeScale;
segments = (segments - startTime) / timeScale;
injections = (injections - startTime) / timeScale;
livetime = (livetime - startTime) / timeScale;
downtimes = (downtimes - startTime) / timeScale;

% determine default time stamp
if isempty(timeStamp),
  timeStamp = [startTimeString ' to ' stopTimeString];
end

% latex compatible channel name
channelName = strrep(channelName, '_', '\_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                sort triggers                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort unclustered triggers
disp('sorting unclustered triggers...');
[ignore, sortedIndices] = sort(triggersRho);
triggersTime = triggersTime(sortedIndices);
triggersFrequency = triggersFrequency(sortedIndices);
triggersQ = triggersQ(sortedIndices);
triggersRho = triggersRho(sortedIndices);

% sort clustered triggers
disp('sorting clustered triggers...');
[ignore, sortedIndices] = sort(clustersRho);
clustersTime = clustersTime(sortedIndices);
clustersFrequency = clustersFrequency(sortedIndices);
clustersQ = clustersQ(sortedIndices);
clustersRho = clustersRho(sortedIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              threshold triggers                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold unclustered triggers
disp('thresholding unclustered triggers...');
keepIndices = find(triggersRho > rho0);
numberOfTriggers = length(keepIndices);
triggersTime = triggersTime(keepIndices);
triggersFrequency = triggersFrequency(keepIndices);
triggersQ = triggersQ(keepIndices);
triggersRho = triggersRho(keepIndices);

% threshold clustered triggers
disp('thresholding clustered triggers...');
keepIndices = find(clustersRho > rho0);
numberOfClusters = length(keepIndices);
clustersTime = clustersTime(keepIndices);
clustersFrequency = clustersFrequency(keepIndices);
clustersQ = clustersQ(keepIndices);
clustersRho = clustersRho(keepIndices);

% identify trigger subsets by snr threshold
triggers10Indices = find(triggersRho > 10);
triggers20Indices = find(triggersRho > 20);
[ignore, triggersLoudestIndex] = max(triggersRho);

% identify cluster subsets by snr threshold
clusters10Indices = find(clustersRho > 10);
clusters20Indices = find(clustersRho > 20);
[ignore, clustersLoudestIndex] = max(clustersRho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              tile signal space                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('tiling signal space...');

% tile signal space
transientFactor = 4;
tiling = wtile(blockDuration, qRange, frequencyRange, sampleFrequency, ...
               maximumMismatch, highPassCutoff, lowPassCutoff, ...
               whiteningDuration, transientFactor);

% construct distribution of tile properties
tilesQ = [];
tilesFrequency = [];
tilesIndependents = [];
for plane = 1 : tiling.numberOfPlanes,
  for row = 1 : tiling.planes{plane}.numberOfRows,
    tilesQ = [tilesQ ...
              tiling.planes{plane}.q];
    tilesFrequency = [tilesFrequency ...
                      tiling.planes{plane}.rows{row}.frequency];
    tilesIndependents = [tilesIndependents ...
                        tiling.planes{plane}.rows{row}.numberOfIndependents];
  end
end
tilesIndependents = 1.5 * tilesIndependents / tiling.duration;
tilesNumberOfBins = length(tilesQ);

% construct coordinate grid for q-frequency distributions
gridQ= log2(unique(sort(tilesQ)));
gridFrequency = log2(unique(sort(tilesFrequency)));
[GridQ, GridFrequency] = meshgrid(gridQ, gridFrequency);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           create output directory                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output directory
unix(['mkdir -p ' outputDirectory]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             setup figure handle                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate plots
generatePlots = true;
% use version 0 printing
printVers = 0;
figureHandle = figure;
%figureHandle = figure('visible', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      unclustered trigger rate vs. time                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger rate vs. time...');

% determine unclustered trigger rate vs. time
binWidth = 60 / timeScale;
binEdges = 0 : binWidth : analysisDuration / timeScale;
binCenters = binWidth / 2 : binWidth : analysisDuration / timeScale - binWidth / 2;
numberOfBins = length(binCenters);
if isempty(livetime),
  binLivetime = zeros(size(binCenters));
else
  binLivetime = ...
      sum(max(0, min(ones(numberOfLivetimes, 1) * binEdges(2 : end), ...
                     livetime(:, 2) * ones(1, numberOfBins)) - ...
                 max(ones(numberOfLivetimes, 1) * binEdges(1 : end - 1), ...
                     livetime(:, 1) * ones(1, numberOfBins))), 1);
end
binLivetime(find(binLivetime == 0)) = 1;
if isempty(triggersTime),
  triggersRate = zeros(length(binCenters), 1);
  triggersRateMinimum = zeros(length(binCenters), 1);
  triggersRateMaximum = zeros(length(binCenters), 1);
else
  triggersCount = histc(triggersTime, binEdges);
  triggersRate = triggersCount(1 : end - 1) ./ (binLivetime * timeScale).';
  triggersRateMinimum = triggersRate - ...
                        sqrt(triggersCount(1 : end - 1)) ./ ...
                        (binLivetime * timeScale).';
  triggersRateMaximum = triggersRate + ...
                        sqrt(triggersCount(1 : end - 1)) ./ ...
                        (binLivetime * timeScale).';
  triggersRateMinimum = max(1e-2, triggersRateMinimum);
  triggersRateMaximum = min(1e+2, triggersRateMaximum);
end

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale 1e-2 1e2])
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles = semilogy([binCenters; binCenters;], ...
                   [triggersRateMinimum.'; triggersRateMaximum.';], ...
                   'k-', ...
                   binCenters, triggersRate.', 'bo');
set(handles(1 : end - 1), 'LineWidth', 1);
set(handles(end), 'MarkerSize', 5, ...
                  'LineWidth', 0.5, ...
                  'MarkerEdgeColor', [0.0 0.0 0.7], ...
                  'MarkerFaceColor', [0.0 0.0 1.0]);
hold off;

% label figure
xlabel(['Time [' timeLabel ']']);
ylabel('Trigger rate [Hz]');
title([channelName ' unclustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_time_rate_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       clustered trigger rate vs. time                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger rate vs. time...');

% determine clustered trigger rate vs. time
binWidth = 60 / timeScale;
binEdges = 0 : binWidth : analysisDuration / timeScale;
binCenters = binWidth / 2 : binWidth : analysisDuration / timeScale - binWidth / 2;
numberOfBins = length(binCenters);
if isempty(livetime),
  binLivetime = zeros(size(binCenters));
else
  binLivetime = ...
      sum(max(0, min(ones(numberOfLivetimes, 1) * binEdges(2 : end), ...
                     livetime(:, 2) * ones(1, numberOfBins)) - ...
                 max(ones(numberOfLivetimes, 1) * binEdges(1 : end - 1), ...
                     livetime(:, 1) * ones(1, numberOfBins))), 1);
end
binLivetime(find(binLivetime == 0)) = 1;
if isempty(clustersTime),
  clustersRate = zeros(length(binCenters), 1);
  clustersRateMinimum = zeros(length(binCenters), 1);
  clustersRateMaximum = zeros(length(binCenters), 1);
else
  clustersCount = histc(clustersTime, binEdges);
  clustersRate = clustersCount(1 : end - 1) ./ (binLivetime * timeScale).';
  clustersRateMinimum = clustersRate - ...
                        sqrt(clustersCount(1 : end - 1)) ./ ...
                        (binLivetime * timeScale).';
  clustersRateMaximum = clustersRate + ...
                        sqrt(clustersCount(1 : end - 1)) ./ ...
                        (binLivetime * timeScale).';
  clustersRateMinimum = max(1e-2, clustersRateMinimum);
  clustersRateMaximum = min(1e+2, clustersRateMaximum);
end

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale 1e-2 1e2])
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles = semilogy([binCenters; binCenters;], ...
                   [clustersRateMinimum.'; clustersRateMaximum.';], ...
                   'k-', ...
                   binCenters, clustersRate.', 'bo');
set(handles(1 : end - 1), 'LineWidth', 1);
set(handles(end), 'MarkerSize', 5, ...
                  'LineWidth', 0.5, ...
                  'MarkerEdgeColor', [0.0 0.0 0.7], ...
                  'MarkerFaceColor', [0.0 0.0 1.0]);
hold off;

% label figure
xlabel(['Time [' timeLabel ']']);
ylabel('Trigger rate [Hz]');
title([channelName ' clustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_time_rate_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    unclustered trigger frequency vs. time                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger frequency vs. time...');


% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale 5 11]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'linear');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles1 = plot(triggersTime, ...
                log2(triggersFrequency), ...
                'bo');
handles2 = plot(triggersTime(triggers10Indices), ...
                log2(triggersFrequency(triggers10Indices)), ...
                'go');
handles3 = plot(triggersTime(triggers20Indices), ...
                log2(triggersFrequency(triggers20Indices)), ...
                'ro');
handles4 = plot(triggersTime(triggersLoudestIndex), ...
                log2(triggersFrequency(triggersLoudestIndex)), ...
                'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
set(gca, 'YTick', 5 : 1 : 11);
set(gca, 'YTickLabel', 2.^(5 : 1 : 11));
xlabel(['Time [' timeLabel ']']);
ylabel('Frequency [Hz]');
title([channelName ' unclustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_time_frequency_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     clustered trigger frequency vs. time                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger frequency vs. time...');

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale 5 11]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'linear');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles1 = plot(clustersTime, ...
                log2(clustersFrequency), ...
                'bo');
handles2 = plot(clustersTime(clusters10Indices), ...
                log2(clustersFrequency(clusters10Indices)), ...
                'go');
handles3 = plot(clustersTime(clusters20Indices), ...
                log2(clustersFrequency(clusters20Indices)), ...
                'ro');
handles4 = plot(clustersTime(clustersLoudestIndex), ...
                log2(clustersFrequency(clustersLoudestIndex)), ...
                'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
set(gca, 'YTick', 5 : 1 : 11);
set(gca, 'YTickLabel', 2.^(5 : 1 : 11));
xlabel(['Time [' timeLabel ']']);
ylabel('Frequency [Hz]');
title([channelName ' clustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure to file
wprintfig(figureHandle,outputDirectory,['clustered_time_frequency_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       unclustered trigger snr vs. time                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger snr vs. time...');

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale rho0 1e3]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles1 = semilogy(triggersTime, ...
                    triggersRho, ...
                    'bo');
handles2 = semilogy(triggersTime(triggers10Indices), ...
                    triggersRho(triggers10Indices), ...
                    'go');
handles3 = semilogy(triggersTime(triggers20Indices), ...
                    triggersRho(triggers20Indices), ...
                    'ro');
handles4 = semilogy(triggersTime(triggersLoudestIndex), ...
                    triggersRho(triggersLoudestIndex), ...
                    'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
xlabel(['Time [' timeLabel ']']);
ylabel('Signal to noise ratio');
title([channelName ' unclustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_time_snr_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        clustered trigger snr vs. time                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger snr vs. time...');

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([0 analysisDuration / timeScale rho0 1e3]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
if numberOfDowntimes > 0,
  handles = fill([downtimes'; fliplr(downtimes)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.8 0.8 0.8]);
  set(handles, 'FaceColor', [0.9 0.9 0.9]);
end
if numberOfInjections > 0,
  handles = fill([injections'; fliplr(injections)';], ...
                 [axisScale(3) axisScale(3) axisScale(4) axisScale(4)]', 'y');
  set(handles, 'LineWidth', 2);
  set(handles, 'EdgeColor', [0.9 0.9 0.0]);
  set(handles, 'FaceColor', [1.0 1.0 0.0]);
end
handles1 = semilogy(clustersTime, ...
                    clustersRho, ...
                    'bo');
handles2 = semilogy(clustersTime(clusters10Indices), ...
                    clustersRho(clusters10Indices), ...
                    'go');
handles3 = semilogy(clustersTime(clusters20Indices), ...
                    clustersRho(clusters20Indices), ...
                    'ro');
handles4 = semilogy(clustersTime(clustersLoudestIndex), ...
                    clustersRho(clustersLoudestIndex), ...
                    'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
xlabel(['Time [' timeLabel ']']);
ylabel('Signal to noise ratio');
title([channelName ' clustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_time_snr_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           mask injection segments                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mask injection segments
disp('masking injection segments...');
triggersInjectionIndices = [];
clustersInjectionIndices = [];
for injectionNumber = 1 : numberOfInjections,
  triggersInjectionIndices = [triggersInjectionIndices; ...
    find((triggersTime >= injections(injectionNumber, 1)) & ...
         (triggersTime < injections(injectionNumber, 2)))];
  clustersInjectionIndices = [clustersInjectionIndices; ...
    find((clustersTime >= injections(injectionNumber, 1)) & ...
         (clustersTime < injections(injectionNumber, 2)))];
end

% mask injection triggers
keepIndices = setdiff(1 : numberOfTriggers, triggersInjectionIndices);
numberOfTriggers = length(keepIndices);
triggersTime = triggersTime(keepIndices);
triggersFrequency = triggersFrequency(keepIndices);
triggersQ = triggersQ(keepIndices);
triggersRho = triggersRho(keepIndices);

% mask injection clusters
keepIndices = setdiff(1 : numberOfClusters, clustersInjectionIndices);
numberOfClusters = length(keepIndices);
clustersTime = clustersTime(keepIndices);
clustersFrequency = clustersFrequency(keepIndices);
clustersQ = clustersQ(keepIndices);
clustersRho = clustersRho(keepIndices);

% identify trigger subsets by snr threshold
triggers10Indices = find(triggersRho > 10);
triggers20Indices = find(triggersRho > 20);
if ismember(triggersLoudestIndex, triggersInjectionIndices),
  triggersLoudestIndex = [];
else
  [ignore, triggersLoudestIndex] = max(triggersRho);
end

% identify cluster subsets by snr threshold
clusters10Indices = find(clustersRho > 10);
clusters20Indices = find(clustersRho > 20);
if ismember(clustersLoudestIndex, clustersInjectionIndices),
  clustersLoudestIndex = [];
else
  [ignore, clustersLoudestIndex] = max(clustersRho);
end

% correct analyzed livetime
analysisLivetime = analysisLivetime - injectionTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    unclustered trigger snr vs. frequency                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger snr vs. frequency...');

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([5 11 rho0 1e3]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
handles1 = semilogy(log2(triggersFrequency), ...
                    triggersRho, ...
                    'bo');
handles2 = semilogy(log2(triggersFrequency(triggers10Indices)), ...
                    triggersRho(triggers10Indices), ...
                    'go');
handles3 = semilogy(log2(triggersFrequency(triggers20Indices)), ...
                    triggersRho(triggers20Indices), ...
                    'ro');
handles4 = semilogy(log2(triggersFrequency(triggersLoudestIndex)), ...
                    triggersRho(triggersLoudestIndex), ...
                    'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 2, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 2, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 2, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
set(gca, 'XTick', 5 : 1 : 11);
set(gca, 'XTickLabel', 2.^(5 : 1 : 11));
xlabel('Frequency [Hz]');
ylabel('Signal to noise ratio');
title([channelName ' unclustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_frequency_snr_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     clustered trigger snr vs. frequency                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger snr vs. frequency...');

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([5 11 rho0 1e3]);
axisScale = axis;
set(gca, 'XScale', 'linear');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
handles1 = semilogy(log2(clustersFrequency), ...
                    clustersRho, ...
                    'bo');
handles2 = semilogy(log2(clustersFrequency(clusters10Indices)), ...
                    clustersRho(clusters10Indices), ...
                    'go');
handles3 = semilogy(log2(clustersFrequency(clusters20Indices)), ...
                    clustersRho(clusters20Indices), ...
                    'ro');
handles4 = semilogy(log2(clustersFrequency(clustersLoudestIndex)), ...
                    clustersRho(clustersLoudestIndex), ...
                    'ko');
set(handles1, 'MarkerSize', 4, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.0 0.7], ...
              'MarkerFaceColor', [0.0 0.0 1.0]);
set(handles2, 'MarkerSize', 6, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.0 0.5 0.0], ...
              'MarkerFaceColor', [0.0 1.0 0.0]);
set(handles3, 'MarkerSize', 8, ...
              'LineWidth', 0.5, ...
              'MarkerEdgeColor', [0.3 0.0 0.0], ...
              'MarkerFaceColor', [1.0 0.0 0.0]);
set(handles4, 'MarkerSize', 20, 'LineWidth', 4);
hold off;

% label figure
set(gca, 'XTick', 5 : 1 : 11);
set(gca, 'XTickLabel', 2.^(5 : 1 : 11));
xlabel('Frequency [Hz]');
ylabel('Signal to noise ratio');
title([channelName ' clustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_frequency_snr_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         unclustered frequency vs. q                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger frequency vs. q...');

% determine q-frequency distribution of unclustered triggers
binCounts = zeros(1, tilesNumberOfBins);
for binNumber = 1 : tilesNumberOfBins,
  binCounts(binNumber) = length(find(...
      (abs(triggersQ ./ tilesQ(binNumber) - 1) < 0.01) & ...
      (abs(triggersFrequency ./ tilesFrequency(binNumber) - 1) < 0.01)));
end
binCounts = binCounts ./ ...
            (analysisLivetime * tilesIndependents * exp(-rho0^2 / 2));
GridCounts = griddata(log2(tilesQ), log2(tilesFrequency), binCounts, ...
                      GridQ, GridFrequency, 'nearest');
warning off;
GridCounts = log2(GridCounts);
warning on;

% plot figure
clf;
set(gca, 'FontSize', 18);
colormapScale = [0 7];
imagesc(gridQ, gridFrequency, GridCounts, colormapScale);
imagePosition = [0.30 0.14 0.60 0.76];
set(gca, 'Position', imagePosition);
set(gca, 'YDir', 'normal');
set(gca, 'XTick', 2 : 7);
set(gca, 'XTickLabel', 2.^(2 : 7));
set(gca, 'YTick', 5 : 11);
set(gca, 'YTickLabel', 2.^(5 : 11));
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.01 0.025]);

% label figure
xlabel('Q');
ylabel('Frequency [Hz]');
title([channelName ' unclustered triggers SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% add colorbar
colorbarPosition = [0.14 0.14 0.02 0.76];
subplot('position', colorbarPosition);
set(gca, 'FontSize', 18);
colorbarmap = linspace(min(colormapScale), max(colormapScale), 100).';
imagesc(1, colorbarmap, colorbarmap, colormapScale);
set(gca, 'YDir', 'normal');
set(gca, 'XTick', [])
set(gca, 'YTick', 0 : 7);
set(gca, 'YTickLabel', 2.^(0 : 7));
set(gca, 'TickDir', 'out')
ylabel('Relative trigger excess');

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_q_frequency_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          clustered frequency vs. q                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger frequency vs. q...');

% determine q-frequency distribution of clustered triggers
binCounts = zeros(1, tilesNumberOfBins);
for binNumber = 1 : tilesNumberOfBins,
  binCounts(binNumber) = length(find(...
      (abs(clustersQ ./ tilesQ(binNumber) - 1) < 0.01) & ...
      (abs(clustersFrequency ./ tilesFrequency(binNumber) - 1) < 0.01)));
end
binCounts = binCounts ./ ...
            (analysisLivetime * tilesIndependents * exp(-rho0^2 / 2));
GridCounts = griddata(log2(tilesQ), log2(tilesFrequency), binCounts, ...
                      GridQ, GridFrequency, 'nearest');
warning off;
GridCounts = log2(GridCounts);
warning on;

% plot figure
clf;
set(gca, 'FontSize', 18);
colormapScale = [0 7];
imagesc(gridQ, gridFrequency, GridCounts, colormapScale);
imagePosition = [0.30 0.14 0.60 0.76];
set(gca, 'Position', imagePosition);
set(gca, 'YDir', 'normal');
set(gca, 'XTick', 2 : 7);
set(gca, 'XTickLabel', 2.^(2 : 7));
set(gca, 'YTick', 5 : 11);
set(gca, 'YTickLabel', 2.^(5 : 11));
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.01 0.025]);

% label figure
xlabel('Q');
ylabel('Frequency [Hz]');
title([channelName ' clustered triggers SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% add colorbar
colorbarPosition = [0.14 0.14 0.02 0.76];
subplot('position', colorbarPosition);
set(gca, 'FontSize', 18);
colorbarmap = linspace(min(colormapScale), max(colormapScale), 100).';
imagesc(1, colorbarmap, colorbarmap, colormapScale);
set(gca, 'YDir', 'normal');
set(gca, 'XTick', [])
set(gca, 'YTick', 0 : 7);
set(gca, 'YTickLabel', 2.^(0 : 7));
set(gca, 'TickDir', 'out')
ylabel('Relative trigger excess');

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_q_frequency_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  unclustered trigger rate vs. snr threshold                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger rate vs. snr threshold...');

% determine unclustered trigger rate vs. snr threshold
triggersRate = (numberOfTriggers : -1 : 1) / analysisLivetime;

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([rho0 1e3 10^floor(log10(1 / analysisDuration)) 1e1]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
handles = loglog(triggersRho, triggersRate, 'bo-');
set(handles, 'MarkerSize', 4);
set(handles, 'LineWidth', 2);
hold off;

% label figure
xlabel('Signal to noise ratio threshold');
ylabel('Trigger rate exceeding threshold [Hz]');
title([channelName ' unclustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_snr_rate_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   clustered trigger rate vs. snr threshold                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting cumulative snr distribution...');

% determine clustered trigger rate vs. snr threshold
clustersRate = (numberOfClusters : -1 : 1) / analysisLivetime;

% set figure axis
clf;
set(gca, 'FontSize', 18);
axis([rho0 1e3 10^floor(log10(1 / analysisDuration)) 1e1]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;

% plot figure
hold on;
handles = loglog(clustersRho, clustersRate, 'bo-');
set(handles, 'MarkerSize', 4);
set(handles, 'LineWidth', 2);
hold off;

% label figure
xlabel('Signal to noise ratio threshold');
ylabel('Trigger rate exceeding threshold [Hz]');
title([channelName ' clustered triggers with SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_snr_rate_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         unclustered autocorrelogram                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting unclustered trigger autocorrelogram...');

% determine unclustered autocorrelogram
timeResolution = 0.01;
frequencyResolution = min(0.002, 100 / analysisDuration);
minimumFrequency = 0.01;
maximumFrequency = 10;
binEdges = (0 : timeResolution : 1 / frequencyResolution) / timeScale;
blockDuration = 1 / frequencyResolution;
fourierLength = 2^nextpow2(blockDuration / timeResolution);
zeroPadLength = fourierLength - length(binEdges);
minimumBlockOverlap = 1 / timeResolution;
if analysisDuration <= blockDuration,
  overlapDuration = 0;
else
  numberOfBlocks = ceil((analysisDuration - minimumBlockOverlap) / ...
                        (blockDuration - minimumBlockOverlap));
  overlapDuration = (analysisDuration - blockDuration * numberOfBlocks) / ...
                    (1 - numberOfBlocks);
end
blockStartTimes = (0 : blockDuration - overlapDuration : ...
                   analysisDuration - blockDuration) / timeScale;
autocorrelogram = zeros(1, fourierLength);
for blockStartTime = blockStartTimes,
  blockStopTime = blockStartTime + (1 / frequencyResolution) / timeScale;
  blockIndices = find((triggersTime >= blockStartTime) & ...
                      (triggersTime < blockStopTime));
  counts = hist(triggersTime(blockIndices) - blockStartTime, binEdges);
  counts = [counts zeros(1, zeroPadLength)];
  autocorrelogram = autocorrelogram + abs(fft(counts)).^2;
end
frequencies = (0 : fourierLength - 1) / (timeResolution * fourierLength);
frequencyResolution = frequencies(2);
autocorrelogram = autocorrelogram * ...
                  (1 / numberOfTriggers) * ...
                  (1 / frequencyResolution) * ...
                  (blockDuration - overlapDuration) / blockDuration;
keepIndices = find((frequencies >= minimumFrequency) & ...
                   (frequencies <= maximumFrequency));
frequencies = frequencies(keepIndices);
autocorrelogram = autocorrelogram(keepIndices);

% plot figure
clf;
set(gca, 'FontSize', 18);
handles = loglog(frequencies, autocorrelogram, 'b-');

% set figure axis
defaultAxis = axis;
axis([minimumFrequency maximumFrequency defaultAxis(3:4)]);
grid on;

% label figure
xlabel('Frequency [Hz]');
ylabel('Number of trigger pairs [Hz^{-1}]');
title([channelName ' unclustered triggers SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['unclustered_autocorrelogram_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          clustered autocorrelogram                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report status
disp('plotting clustered trigger autocorrelogram...');

% determine clustered autocorrelogram
timeResolution = 0.01;
frequencyResolution = min(0.002, 100 / analysisDuration);
minimumFrequency = 0.01;
maximumFrequency = 10;
binEdges = (0 : timeResolution : 1 / frequencyResolution) / timeScale;
blockDuration = 1 / frequencyResolution;
fourierLength = 2^nextpow2(blockDuration / timeResolution);
zeroPadLength = fourierLength - length(binEdges);
minimumBlockOverlap = 1 / timeResolution;
if analysisDuration <= blockDuration,
  overlapDuration = 0;
else
  numberOfBlocks = ceil((analysisDuration - minimumBlockOverlap) / ...
                        (blockDuration - minimumBlockOverlap));
  overlapDuration = (analysisDuration - blockDuration * numberOfBlocks) / ...
                    (1 - numberOfBlocks);
end
blockStartTimes = (0 : blockDuration - overlapDuration : ...
                   analysisDuration - blockDuration) / timeScale;
autocorrelogram = zeros(1, fourierLength);
for blockStartTime = blockStartTimes,
  blockStopTime = blockStartTime + (1 / frequencyResolution) / timeScale;
  blockIndices = find((clustersTime >= blockStartTime) & ...
                      (clustersTime < blockStopTime));
  counts = hist(clustersTime(blockIndices) - blockStartTime, binEdges);
  counts = [counts zeros(1, zeroPadLength)];
  autocorrelogram = autocorrelogram + abs(fft(counts)).^2;
end
frequencies = (0 : fourierLength - 1) / (timeResolution * fourierLength);
frequencyResolution = frequencies(2);
autocorrelogram = autocorrelogram * ...
                  (1 / numberOfClusters) * ...
                  (1 / frequencyResolution) * ...
                  (blockDuration - overlapDuration) / blockDuration;
keepIndices = find((frequencies >= minimumFrequency) & ...
                   (frequencies <= maximumFrequency));
frequencies = frequencies(keepIndices);
autocorrelogram = autocorrelogram(keepIndices);

% plot figure
clf;
set(gca, 'FontSize', 18);
handles = loglog(frequencies, autocorrelogram, 'b-');

% set figure axis
defaultAxis = axis;
axis([minimumFrequency maximumFrequency defaultAxis(3:4)]);
grid on;

% label figure
xlabel('Frequency [Hz]');
ylabel('Number of trigger pairs [Hz^{-1}]');
title([channelName ' clustered triggers SNR > ' num2str(rho0)]);
handle = axes('Position', get(gca, 'Position'));
set(handle, 'FontSize', 18, 'YTick', [], 'XTick', [], ...
            'Color', 'none', 'YAxisLocation', 'right');
ylabel(timeStamp);

% print figure
wprintfig(figureHandle,outputDirectory,['clustered_autocorrelogram_' num2str(rho0)],generatePlots,printVers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          close all figures and quit                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all figures
close all;

% return to calling function
return;
