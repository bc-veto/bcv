function candidate = wfindcandidate(clusters,triggersClustered,time,debugLevel)
% WFINDCANDIDATE find a candidate event from triggers and return the entire
% cluster in which that trigger resides.
%
% Authors:
% Jameson Rollins <jrollins@phys.columbia.edu>

if (nargin < 2) || isempty(time),
  time = [];
end

if isempty(time),
  wlog(debugLevel, 1, '    finding loudest trigger\n');
  [channel, index, value] = findloudest(clusters);
else
  wlog(debugLevel, 1, '    finding trigger closest to %f\n', time);
  [channel, index, value] = findtime(clusters,time);
end

if isempty(channel),
  error('no candidate channel found??');
end

% get cluster number in which candidate trigger resides
clusterNumber = clusters{channel}.clusterNumber(index);

% get all cluster triggers
candidate = wpullcluster(triggersClustered,clusters,channel,clusterNumber);

% report info on candidate candidate
wlog(debugLevel, 1, '    candidate cluster:\n');
wlog(debugLevel, 1, '      %-26s %s\n', 'channel:', candidate.channelName);
wlog(debugLevel, 1, '      %-26s %f\n', 'time:', candidate.clusterTime);
wlog(debugLevel, 1, '      %-26s %f\n', 'frequency:', candidate.clusterFrequency);
wlog(debugLevel, 1, '      %-26s %f\n', 'duration:', candidate.clusterDuration);
wlog(debugLevel, 1, '      %-26s %f\n', 'bandwidth:', candidate.clusterBandwidth);
wlog(debugLevel, 1, '      %-26s %d\n', 'triggers:', candidate.clusterSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [channel,index,value] = findtime(clusters,time)
% find cluster closest to time

  channel = [];
  value = [];
  index = [];
  P = Inf;
  for channelNumber = 1 : length(clusters),
    % trigger time closest to time
    [M,I] = min(abs(time - clusters{channelNumber}.time));
    if (M < P),
      P = M;
      channel = channelNumber;
      index = I;
      value = clusters{channelNumber}.time(I);
    end
  end

  return


function [channel,index,value] = findloudest(clusters)
% find the overall loudest cluster

  % FIXME: this should be trigger with loudest event in second loudest channel
  channel = [];
  value = 0;
  index = [];
  for channelNumber = 1 : length(clusters),
    % loudest cluster
    [M,I] = max(clusters{channelNumber}.normalizedEnergy);
    if (M > value),
      channel = channelNumber;
      index = I;
      value = M;
    end
  end

  return
