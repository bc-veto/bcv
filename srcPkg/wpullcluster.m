function cluster = wpullcluster(triggers,clusters,channel,clusterNumber)
% WPULLCLUSTER return triggers from a single cluster
%
% Authors:
% Jameson Rollins <jrollins@phys.columbia.edu>

% get the trigger cluster indices
triggerIndices = find(triggers{channel}.clusterNumber == clusterNumber);

% copy all the triggers from the cluster
cluster = wcopyevents(triggers{channel},triggerIndices);

clusterIndex = find(clusters{channel}.clusterNumber == clusterNumber);

cluster.channelNumber = channel;

% add cluster info to the event structure
cluster.clusterNumber = clusterNumber;
cluster.clusterSize = clusters{channel}.size(clusterIndex);
cluster.clusterTime = clusters{channel}.time(clusterIndex);
cluster.clusterFrequency = clusters{channel}.frequency(clusterIndex);
cluster.clusterDuration = clusters{channel}.duration(clusterIndex);
cluster.clusterBandwidth = clusters{channel}.bandwidth(clusterIndex);
cluster.clusterNormalizedEnergy = clusters{channel}.normalizedEnergy(clusterIndex);
