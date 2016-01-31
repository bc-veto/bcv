function clusters = wcoincidecluster(events, clusters, debugLevel)
% WCOINCIDECLUSTER weed out clusters not in trigger list
%
% WCOINCIDECLUSTER(events, clusters, debugLevel)
%
% See also WCOINCIDE, WCLUSTER

% Jameson Rollins <jrollins@phys.columbia.edu>

% $Id: wcoincidecluster.m 1544 2009-03-16 18:57:15Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 3, nargin));

% apply default arguments
if (nargin < 3) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
events = wmat2cell(events);
clusters = wmat2cell(clusters);

% force one dimensional cell arrays
events = events(:);
clusters = clusters(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(events);

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(events{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('events argument is not a discrete Q transform event structure');
  end
  if ~strcmp(clusters{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('clusters argument is not a discrete Q transform event structure');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       check for triggers in clusters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for channelNumber = 1 : numberOfChannels,

  % determine the coincident cluster indices from the unique event cluster
  % numbers
  coincidentIndices = unique(events{channelNumber}.clusterNumber);

  clusters{channelNumber} = wcopyevents(clusters{channelNumber}, ...
                                       coincidentIndices);
  
  % write the original cluster numbers to a clusterNumber element
  clusters{channelNumber}.clusterNumber = coincidentIndices;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return consistent clusters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
