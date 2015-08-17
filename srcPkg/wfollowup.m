function [event, skymap] = ...
        wfollowup(data, coefficients, tiling, blockStartTime, coordinate, ...
                  candidateCluster, channelNames, parameters, debugLevel)
% WFOLLOWUP Omega Bayesian/coherent followup code
%
% WFOLLOWUP handles execution of all follow-up analyses on coincident single
% detector triggers.
%
% usage: [event, skymap] = ...
%        wfollowup(data, coefficients, tiling, blockTime, coordinate, ...
%                  candidateCluster, channelNames, parameters, debugLevel)
%
%   data                  cell array of input time series
%   coefficients          cell array of frequency domain filter coefficients
%   tiling                q tiling structure
%   blockStartTime        block start time
%   coordinate            geocentric sky position coordinate
%   candidateCluter       structure of cluster triggers to search over
%   channelNames          channel names
%   parameters            parameter structure

% Authors:
% Jameson Rollins <jrollins@phys.columbia.edu>
% Antony Searle <antony.searle@anu.edu.au>

% $Id: wfollowup.m 2336 2009-09-25 20:57:23Z acsearle $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             analysis parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(9, 9, nargin));

% default Bayesian parameters
durationInflation = 0.030; % ms

% extract needed parameters from parameters structure
maxFollowTriggers = parameters.maxFollowTriggers;
xCoherentCheck = parameters.xCoherentCheck;

% number of channels
numberOfChannels = length(channelNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      create the output event structure                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event.id = 'Discrete Q-transform event structure';
event.time = candidateCluster.clusterTime;
event.frequency = candidateCluster.clusterFrequency;
event.duration = candidateCluster.clusterDuration;
event.bandwidth = candidateCluster.clusterBandwidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         determine sky directions                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if coordinate specified, get theta and phi from coordinate
if ~isempty(coordinate),

  % otherwise get theta and phi from coordinate
  theta = coordinate(1);
  phi = coordinate(2);

else

  % if no coordinate is specified, get coordinates from projection
  [theta, phi] = wsinusoidalprojection(180*2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              bayesian analysis                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wlog(debugLevel, 1, '  computing bayesian posteriors\n');

% construct the input data to wposteriors
% FIXME: comment the code in this section
% FIXME: get a per-detector time estimate
% FIXME: is padding with durationInflation right?

logSkymap = zeros(1, length(theta));
logGlitches = zeros(1, numberOfChannels);

% determine trigger energy ranking (triggers should be already sorted, but just
% in case...)
[dummy, energySortIndices] = ...
    sort(candidateCluster.normalizedEnergy(:),1,'descend');
clear dummy;

% loop over the triggers (limited to max specified)
for index = energySortIndices(1:min(candidateCluster.clusterSize,maxFollowTriggers))',

    time = candidateCluster.time(index);
    frequency = candidateCluster.frequency(index);
    duration = candidateCluster.duration(index);
    bandwidth = candidateCluster.bandwidth(index);
    q = candidateCluster.q(index);

    wlog(debugLevel, 2, '    cluster trigger %d\n', index);
    wlog(debugLevel, 2, '      %-26s %f\n', 'time:', time);
    wlog(debugLevel, 2, '      %-26s %f\n', 'frequency:', frequency);
    wlog(debugLevel, 2, '      %-26s %f\n', 'duration:', duration);
    wlog(debugLevel, 2, '      %-26s %f\n', 'bandwidth:', bandwidth);
    wlog(debugLevel, 2, '      %-26s %f\n', 'q:', q);

    minFrequency = frequency - bandwidth / 2;
    maxFrequency = frequency + bandwidth / 2;
    minTime = time - duration / 2 - blockStartTime - durationInflation;
    maxTime = time + duration / 2 - blockStartTime + durationInflation;

    wlog(debugLevel, 3, '        [ %g, %g ] s\n', minTime, maxTime);
    wlog(debugLevel, 3, '        [ %g, %g ] Hz\n', minFrequency, maxFrequency);

    qIndex = find(tiling.qs == q);
    rowIndex = find(tiling.planes{qIndex}.frequencies == frequency);

    times = [];

    sampleFrequency = 16384; % crank the rate up

    for channelNumber = 1:numberOfChannels,

	% deduce size of FFT buffer from one-sided data FFT
        xSw{channelNumber} = zeros(1, (length(data{channelNumber}) - 1) * 2 * sampleFrequency / tiling.sampleFrequency);

        xSw{channelNumber}(tiling.planes{qIndex}.rows{rowIndex}.dataIndices) = ...
            tiling.planes{qIndex}.rows{rowIndex}.window .* ...
            data{channelNumber}(tiling.planes{qIndex}.rows{rowIndex}.dataIndices);

        xSw{channelNumber} = ifft(xSw{channelNumber}) * length(xSw{channelNumber}) * sqrt(2);

        wSw(channelNumber) = ...
            sum((tiling.planes{qIndex}.rows{rowIndex}.window .* ...
             coefficients{channelNumber}(tiling.planes{qIndex}.rows{rowIndex}.dataIndices)) .^2);

        % compute the normalized matched filter

        xSwqwSw{channelNumber} = xSw{channelNumber} ./ sqrt(wSw(channelNumber));

	if debugLevel > 2
            % check the normalization; these should be approxeq 1
            realmean = mean(real(xSwqwSw{channelNumber}).^2);
            imagmean = mean(imag(xSwqwSw{channelNumber}).^2);
            wlog(debugLevel, 3, '        channel %d normalization: %g, %g\n', channelNumber, realmean, imagmean);
	end

	% use the normalized matched filter to make a 95% confidence region

        minTimeI = floor(minTime * sampleFrequency);
        maxTimeI = ceil(maxTime * sampleFrequency);

        y = 0.25 * abs(xSwqwSw{channelNumber}(minTimeI:maxTimeI)).^2;
        [dummy, k] = sort(y, 'descend');
        z = cumsum(exp(y(k) - max(y)));
        z = z / z(end);
        z = find(z > 0.99, 1);
        z = y >= y(k(z));
        minT = (find(z, 1) -1 - 1 + minTimeI) / sampleFrequency;
        maxT = (find(z, 1, 'last') -1 + 1 + minTimeI) / sampleFrequency;
    
        times = [ times, minT, maxT ];

    end

    dTau = 0.1 * duration;

    wlog(debugLevel, 3, '        %-25s%g, %g\n', 'dTau (s, samples):', dTau, dTau * sampleFrequency);

    times = [ times, minTime, maxTime, dTau ];

    [logSkymapI, logGlitchesI] = wposteriorsgs(channelNames, sampleFrequency, wSw, xSw, [theta' ; phi'], times);    
    
    % uncomment to interpret the event as any combination of signals and glitches
    %logSkymapI = logsumexp(logSkymapI, sum(logGlitchesI)) - log(2);

    % don't delete me, i'm useful!
    %figure(index);
    %whybridskyplot([theta, phi, logSkymapI']);
    
    logSkymap = logSkymap + logSkymapI;
    logGlitches = logGlitches + logGlitchesI;

    wlog(debugLevel, 2, '      %-26s %f\n', 'logSignal:', log(mean(exp(logSkymapI - max(logSkymapI)))) + max(logSkymapI));
    wlog(debugLevel, 2, '      %-26s %f\n', 'logGlitch:', sum(logGlitchesI));

end

% don't delete me, i'm helpful
%figure;
%whybridskyplot([theta, phi, logSkymap']);
%figure;
%whybridskyplot([theta, phi, exp(logSkymap' - max(logSkymap'))]);

% unpack posteriors into constituent pieces
event.logSignal = log(mean(exp(logSkymap - max(logSkymap)))) + max(logSkymap);
event.logGlitch = sum(logGlitches);
[dummy, i] = max(logSkymap);
event.modeTheta = theta(i);
event.modePhi = phi(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           create bayesian skymap                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the ranking of the pixels from highest to lowest
[logSkymap, index] = sort(logSkymap(:),1,'descend');
theta = theta(index);
phi = phi(index);

% normalize the skymap
%logSkymap = logSkymap - logSkymap(1);
%skymap = exp(logSkymap);
%skymap = skymap / sum(skymap(:));
%skymap = skymap(:);

skymap = [theta, phi, logSkymap];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          compute event statistics                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the probabilities for model selection.
% probSignal is "the" detection statistic, ranging from 0 (sure not a
% signal) to 1/3 (data has not informed choice between three hypotheses) to
% 0.95 (95% confident a signal is present) to 1 (sure that a signal is 
% present)
event.probSignal = ...
    1 / (1 + exp(event.logGlitch - event.logSignal) + exp(-event.logSignal));
event.probGlitch = ...
    1 / (exp(event.logSignal - event.logGlitch) + 1 + exp(-event.logGlitch));
event.probNoise  = ...
    1 / (exp(event.logSignal) + exp(event.logGlitch) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           report event properties                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% report Bayesian event properties
wlog(debugLevel, 1, '    event information:\n');
wlog(debugLevel, 1, '      %-26s %f\n', 'modeTheta:', event.modeTheta);
wlog(debugLevel, 1, '      %-26s %f\n', 'modePhi:', event.modePhi);
wlog(debugLevel, 1, '      %-26s %f\n', 'probGlitch:', event.probGlitch);
wlog(debugLevel, 1, '      %-26s %f\n', 'probSignal:', event.probSignal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
