function [significants, clusters] = ...
    wcluster(significants, clusterMethod, ...
             clusterParameter1, clusterParameter2, clusterParameter3, ...
             distanceMetric, durationInflation, bandwidthInflation, ...
             debugLevel)
% WCLUSTER clustering of significant Q transform tiles
%
% WCLUSTER identifies clusters of related significant Q transform tiles in order
% to improve the detection of gravitational-wave bursts that are extended in
% time and/or frequency.
%
% WCLUSTER permits two choices of clustering methods: density and hierarchical.
% The input arguments depends upon which method is requested.  By default, the
% 'hierarchical' method is assumed.
%
% usage:
%
%        [significants, clusters] = ...
%            wcluster(significants, 'density', ...
%                     clusterRadius, clusterDensity, clusterSingles, ...
%                     distanceMetric, durationInflation, bandwidthInflation, ...
%                     debugLevel);
% 
%        [significants, clusters] = ...
%            wcluster(significants, 'hierarchical', ...
%                     clusterLinkage, clusterCriterion, clusterThreshold, ...
%                     distanceMetric, durationInflation, bandwidthInflation, ...
%                     debugLevel);
% 
% The following input arguments are independent of clustering method.
%
%   significants         cell array of significant tile properties
%   clusterMethod        clustering method 'hiearichal' or 'density'
%   distanceMetric       choice of metric for computing distance
%   durationInflation    multiplicative scale factor for duration
%   bandwidthInflation   multiplicative scale factor for bandwidth
%   debugLevel           verboseness of debug output
%
% The following input arguments are specific to the density method.
%
%   clusterRadius        length scale for clustering
%   clusterDensity       required number of tiles in clusterRadius
%   clusterSingles       
%
% The following input arguments are specific to the hierarchical method.
%
%   clusterLinkage       linkage method for merging clusters
%   clusterCriterion     criteria for identifying clusters
%   clusterThreshold     threshold for identifying clusters
%
% The following output arguments are returned.
%
%   significants        cell array of clustered significant tile properties
%   clusters            cell array of aggregate cluster properties
%
% For both density and hierarchical clustering, distances between all pairs of
% tiles are first computed using the WDISTANCE function and the specified
% distance metric and inflation parameters.
%
% The distance metric specified how distances are computed between tiles, and
% may be any of the methods understood by the WDISTANCE function.  By default,
% 'integratedMismatch' is assumed.
%
%   'pointMismatch'        second order approximation to mismatch function
%                          evaluated at center point between tiles
%   'integratedMismatch'   second order approximation to mismatch function
%                          integrated between tiles
%   'logMismatch'          negative log of the overlap between two tiles
%   'euclidean'            normalized time frequency distance between tiles
%   'modifiedEuclidean'    normalized time frequency distance between tiles
%                          modified to give greater weight to frequency
%
% Density clustering:
%
%   Density based clustering first identifies the number of tiles within a
%   specified clustering radius of a given tile.  If the number of tiles within
%   this radius meets a threshold cluster density, then all of the tiles within
%   this radius become part of the same cluster.  Clusters are constructed by
%   recursively applying this same test to all of the other tiles witin the
%   clustering radius, until a complete set of clusters have been formed.  By
%   default, both clusterRadius and clusterDensity are set to unity.
% 
%   Density based clustering does not in general assign all tiles to clusters,
%   but if the cluster singles argument is non zero, single tiles are included
%   in the output cluster list.  Otherwise, single tiles are not included, and
%   zeros are reported for their cluster properties in the output significants
%   structure.  By deafult, clusterSingles is zero.
%
% Hierarchical clustering:
%
%   In hierarchical clustering, the distances returned by WDISTANCE are used by
%   the LINKAGE function to construct a hierarchical cluster tree by
%   successively linking together the next closest pair of tiles or existing
%   clusters.  Finally, clusters are identified by cutting the hieararchical
%   cluster tree based on the specified cretirion and threshold.
%
%   The cluster linkage method defines how intercluster distance is computed,
%   and may be any method understood by the LINKAGE function.  By default,
%   'single' linkage is assumed.
%
%     'single'               nearest distance between cluster members
%     'complete'             furthest distance between cluster members
%     'average'              average distance between cluster members
%     'centroid'             center of mass distance between clusters
%     'ward'                 inner squared distance between clusters
%
%   The cluster criterion specifies how to cut the hierarchical cluster tree to
%   form clusters, and may be either of the two methods understood by the
%   CLUSTER function.  By deafult, 'inconsistent' is assumed.
%
%     'inconsistent'         threshold on inconsistency between of cluster
%                            distance and average subcluster difference
%     'distance'             threshold on absolute distance between clusters
%
%   The clusterThreshold corresponds to the specified clusterCriterion, and has
%   a default threshold of unity for both criteria.
%
% The specified duration and bandwidth inflation factors are used by some of
% distance metrics to compute the distance between tiles.  They are also used
% when computing the total energies of clusters.  By default, they are set to
% unity.
%
% The optional debugLevel argument provides control over the verboseness of
% diagnostic output.  By default, debugLevel is set to unity.
%
% The returned significant tile structure contains the following additional
% fields.
%
%   clusterNumber            cluster identification number
%   clusterSize              number of tiles in cluster
%   clusterNormalizedEnergy  total normalized energy of cluster
%
% For coherent transform data, the following field is also returned
%
%   clusterIncoherentEnergy  total incoherent energy of cluster
%
% The returned cell array of cluster structures contains the following
% fields.
%
%   size              number of tiles in cluster
%   time              characteristic center time of cluster
%   frequency         characteristic center frequency of cluster
%   duration          characteristic duration of cluster
%   bandwidth         characteristic bandwidth of cluster
%   normalizedEnergy  total normalized energy of cluster
%
% For coherent transform data, the following field is also returned
%
%   incoherentEnergy  total incoherent energy of cluster
%
% See also WTILE, WTRANSFORM, WTHRESHOLD, WSELECT, WDISTANCE, LINKAGE,
% INCONSISTENT, and CLUSTER.

% Rubab Khan <rmk2109@columbia.edu>
% Shourov Chatterji <shourov@ligo.caltech.edu>

% $Id: wcluster.m 1644 2009-04-01 16:25:09Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            hard coded parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum allowed number of recursions for density based clustering
maximumRecursions = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 9, nargin));

% apply default arguments
if (nargin < 2) || isempty(clusterMethod),
  clusterMethod = 'hierarchical';
end
switch clusterMethod,
  case 'density',
    if (nargin < 3) || isempty(clusterParameter1),
      clusterRadius = 1.0;
    else
      clusterRadius = clusterParameter1;
    end
    if (nargin < 4) || isempty(clusterParameter2),
      clusterDensity = 1.0;
    else
      clusterDensity = clusterParameter2;
    end
    if (nargin < 5) || isempty(clusterParameter3),
      clusterSingles = 0;
    else
      clusterSingles = clusterParameter3;
    end
  case 'hierarchical',
    if (nargin < 3) || isempty(clusterParameter1),
      clusterLinkage = 'single';
    else
      clusterLinkage = clusterParameter1;
    end
    if (nargin < 4) || isempty(clusterParameter2),
      clusterCriterion = 'inconsistent';
    else
      clusterCriterion = clusterParameter2;
    end
    if (nargin < 5) || isempty(clusterParameter3),
      clusterThreshold = 1.0;
    else
      clusterThreshold = clusterParameter3;
    end
  otherwise,
    error('unknown clustering method %s\n', clusterMethod);
end
if (nargin < 6) || isempty(distanceMetric),
  distanceMetric = 'integratedMismatch';
end
if (nargin < 7) || isempty(durationInflation),
  durationInflation = 1.0;
end
if (nargin < 8) || isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
end
if (nargin < 9) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
significants = wmat2cell(significants);

% force one dimensional cell arrays
significants = significants(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(significants);

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(significants{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              compute distances                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute pairwise distance between tiles
distances = wdistance(significants, distanceMetric, ...
                      durationInflation, bandwidthInflation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      initialize clusters structure                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of cluster structures
clusters = cell(size(significants));

% begin loop over channels
for channelNumber = 1 : numberOfChannels

  % insert structure identification string
  clusters{channelNumber}.id = 'Discrete Q-transform event structure';
  clusters{channelNumber}.channelName = significants{channelNumber}.channelName;

  % output channel name
  significants{channelNumber}.channelName = ...
      regexprep(significants{channelNumber}.channelName, ':.*$', ':CLUSTERED');

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  % number of tiles
  numberOfTiles = length(significants{channelNumber}.time);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                    begin switch on clustering method                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % begin switch on clustering method
  switch clusterMethod,
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          density clustering                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % handle case of density clustering
    case 'density',

      % convert distance structure to matrix
      distanceMatrix = squareform(distances{channelNumber}.distance);

      % handle single tile case
      if numberOfTiles == 1,
        distanceMatrix = 0;
      end

      % initialize the tiles structure
      clear tiles;
      tiles(1 : numberOfTiles) = ...
          struct('clusterNumber', NaN, 'neighborTileNumbers', []);

      % begin loop over tiles
      for tileNumber = 1 : numberOfTiles,

        % find indices of neighboring tiles
        tiles(tileNumber).neighborTileNumbers = ...
            find((distanceMatrix(:, tileNumber) <= clusterRadius) & ...
                 (distanceMatrix(:, tileNumber) ~= 0));

        % if number of tiles within cluster radius exceeds critical density,
        if length(tiles(tileNumber).neighborTileNumbers) + 1 >= clusterDensity,

          % identify tile as potential seed
          tiles(tileNumber).clusterNumber = 0;

          % sort neighboring tiles by increasing normalized energy
          [ignore, sortedNeighborIndices] = ...
              sort(significants{channelNumber}.normalizedEnergy( ...
                  tiles(tileNumber).neighborTileNumbers));
          clear ignore;
          sortedNeighborIndices = fliplr(sortedNeighborIndices);
          tiles(tileNumber).neighborTileNumbers = ...
              tiles(tileNumber).neighborTileNumbers(sortedNeighborIndices);

        % end test for critical density
        end

      % end loop over tiles
      end

      % free distance matrix memory
      clear distanceMatrix;

      % sort significant tiles by increasing normalized energy
      [ignore, sortedTileIndices] = ...
          sort(significants{channelNumber}.normalizedEnergy);
      clear ignore;
      sortedTileIndices = fliplr(sortedTileIndices);

      % initialize cluster number counter
      clusterNumber = 0;

      % begin loop over sorted tiles
      for sortedTileNumber = 1 : numberOfTiles,

        % find tile number corersponding to sorted tile number
        tileNumber = sortedTileIndices(sortedTileNumber);

        % if current tile has not been processed,
        if (tiles(tileNumber).clusterNumber == 0),

          % create a new cluster
          clusterNumber = clusterNumber + 1;

          % assign current tile to new cluster
          tiles(tileNumber).clusterNumber = clusterNumber;

          % set recursion limit to avoid unwanted errors
          previousRecursionLimit = get(0, 'RecursionLimit');
          set(0, 'RecursionLimit', maximumRecursions + 10);

          % find other tiles in the cluster
          tiles = recurse(tiles, tileNumber, maximumRecursions);

          % reset recursion limit
          set(0, 'RecursionLimit', previousRecursionLimit);

          % check for merge with existing cluster
          if tiles(tileNumber).clusterNumber ~= clusterNumber,
            clusterNumber = clusterNumber - 1;
          end
          
        % otherwise, skip to the next tile
        end

      % end loop over sorted tiles
      end

      % collect cluster numbers
      significants{channelNumber}.clusterNumber = [tiles(:).clusterNumber];

      % construct single tile clusters if requested
      singlesIndices = find(isnan(significants{channelNumber}.clusterNumber));
      if clusterSingles,
        significants{channelNumber}.clusterNumber(singlesIndices) = ...
            clusterNumber + (1 : length(singlesIndices));
      else
        significants{channelNumber}.clusterNumber(singlesIndices) = 0;
      end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       hierarchical clustering                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % handle case of hierarchical clustering
    case 'hierarchical',

      % switch on number of significant tiles
      switch numberOfTiles,

        % handle case of zero significant tiles
        case 0,

          % provide empty cluster number list
          significants{channelNumber}.clusterNumber = [];

        % handle case of one significant tile
        case 1,

          % assign single tile to a single cluster
          significants{channelNumber}.clusterNumber = 1;

        % handle case of more than one significant tile
        otherwise,

          % produce heirarchical cluster tree
          tree = linkage(distances{channelNumber}.distance, clusterLinkage);

          % construct clusters from tree
          significants{channelNumber}.clusterNumber = ...
              cluster(tree, 'cutoff', clusterThreshold, ...
                      'criterion', clusterCriterion)';

      % end switch on number of significant tiles
      end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      unknown clustering method                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % handle case of unknown clustering method
    otherwise,
    
      % report error
      error('unknown clustering method %s\n', clusterMethod);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                     end switch on clustering method                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % end switch on clustering method
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         cluster properties                                 %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of clusters
  numberOfClusters = max(significants{channelNumber}.clusterNumber);

  % initialize significant event structure fields
  significants{channelNumber}.clusterSize = ...
      zeros(1, numberOfTiles);
  significants{channelNumber}.clusterNormalizedEnergy = ...
      zeros(1, numberOfTiles);
  if isfield(significants{channelNumber}, 'incoherentEnergy'),
    significants{channelNumber}.clusterIncoherentEnergy = ...
        zeros(1, numberOfTiles);
  end

  % initialize clustered event structure fields
  clusters{channelNumber}.size = zeros(1, numberOfClusters);
  clusters{channelNumber}.time = zeros(1, numberOfClusters);
  clusters{channelNumber}.frequency = zeros(1, numberOfClusters);
  clusters{channelNumber}.duration = zeros(1, numberOfClusters);
  clusters{channelNumber}.bandwidth = zeros(1, numberOfClusters);
  clusters{channelNumber}.normalizedEnergy = zeros(1, numberOfClusters);
  if isfield(significants{channelNumber}, 'incoherentEnergy'),
    clusters{channelNumber}.incoherentEnergy = zeros(1, numberOfClusters);
  end
  
  % begin loop over clusters
  for clusterNumber = 1 : numberOfClusters,

    % find tiles in cluster
    indices = ...
        find(significants{channelNumber}.clusterNumber == clusterNumber);

    % cluster size
    clusters{channelNumber}.size(clusterNumber) = ...
        length(indices);
    significants{channelNumber}.clusterSize(indices) = ...
        clusters{channelNumber}.size(clusterNumber);
    
    % cluster normalized energy
    clusters{channelNumber}.normalizedEnergy(clusterNumber) = ...
        sum(significants{channelNumber}.normalizedEnergy(indices)) * ...
        durationInflation * bandwidthInflation;
    significants{channelNumber}.clusterNormalizedEnergy(indices) = ...
        clusters{channelNumber}.normalizedEnergy(clusterNumber);
    
    % cluster incoherent energy
    if isfield(significants{channelNumber}, 'incoherentEnergy'),
      clusters{channelNumber}.incoherentEnergy(clusterNumber) = ...
          sum(significants{channelNumber}.incoherentEnergy(indices)) * ...
          durationInflation * bandwidthInflation;
        significants{channelNumber}.clusterIncoherentEnergy(indices) = ...
        clusters{channelNumber}.incoherentEnergy(clusterNumber);
    end

    % cluster time
    clusters{channelNumber}.time(clusterNumber) = ...
        sum(significants{channelNumber}.time(indices) .* ...
            significants{channelNumber}.normalizedEnergy(indices)) * ...
        durationInflation * bandwidthInflation ./ ...
        clusters{channelNumber}.normalizedEnergy(clusterNumber);
    
    % cluster frequency
    clusters{channelNumber}.frequency(clusterNumber) = ...
        sum(significants{channelNumber}.frequency(indices) .* ...
            significants{channelNumber}.normalizedEnergy(indices)) * ...
        durationInflation * bandwidthInflation ./ ...
        clusters{channelNumber}.normalizedEnergy(clusterNumber);
    
    % cluster duration
    clusters{channelNumber}.duration(clusterNumber) = ...
        sqrt(sum((significants{channelNumber}.duration(indices).^2 + ...
                  (significants{channelNumber}.time(indices) - ...
                   clusters{channelNumber}.time(clusterNumber)).^2) .* ...
                 significants{channelNumber}.normalizedEnergy(indices)) * ...
             durationInflation * bandwidthInflation ./ ...
             clusters{channelNumber}.normalizedEnergy(clusterNumber));
    
    % cluster bandwidth
    clusters{channelNumber}.bandwidth(clusterNumber) = ...
        sqrt(sum((significants{channelNumber}.bandwidth(indices).^2 + ...
                  (significants{channelNumber}.frequency(indices) - ...
                   clusters{channelNumber}.frequency(clusterNumber)).^2) .* ...
                 significants{channelNumber}.normalizedEnergy(indices)) * ...
             durationInflation * bandwidthInflation ./ ...
             clusters{channelNumber}.normalizedEnergy(clusterNumber));
    
  % end loop over clusters
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               return clusters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            recursion subfunction                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tiles = recurse(tiles, tileNumber, maximumRecursions, recursionNumber)
% RECURSE Auxiliary function for density based clustering algorithm
%
% RECURSE is an auxiliary function for the WCLUSTER density based clustering
% algorithm.  It implements recursive clustering of tiles and should only be
% called by WCLUSTER.
%
% usage: tiles = recurse(tiles, tileNumber, maximumRecursions, recursionNumber);
%
%   tiles               input tile neighbor and cluster structure
%   tileNumber          current tile under test
%   maximumRecursions   limit on allowed recursion depth
%   recursionNumber     current recursion depth
%
%   tiles               updated tile neighbor and cluster structure
%
% See also WCLUSTER.

% verify correct number of input arguments
error(nargchk(2, 4, nargin));

% apply default arguments
if nargin < 3 || isempty(maximumRecursions),
  maximumRecursions = 100;
end
if nargin < 4 || isempty(recursionNumber),
  recursionNumber = 1;
end

% return if recursion limit is exceeded
if (recursionNumber >= maximumRecursions),
  tiles(tileNumber).clusterNumber = 0;
  return;
end

% increment recursion number counter
recursionNumber = recursionNumber + 1;

% determine number of significant tiles
numberOfTiles = length(tiles);

% begin loop over neighboring tiles
for neighborTileNumber = tiles(tileNumber).neighborTileNumbers',

  % if neighbor tile has less than critical density,
  if isnan(tiles(neighborTileNumber).clusterNumber),

    % assign neighbor tile as border tile of current cluster
    tiles(neighborTileNumber).clusterNumber = tiles(tileNumber).clusterNumber;

  % if neighbor tile has critical density and have not been processed,
  elseif tiles(neighborTileNumber).clusterNumber == 0,

    % assign neighbor tile to current cluster
    tiles(neighborTileNumber).clusterNumber = tiles(tileNumber).clusterNumber;

    % continue to recursively build the cluster
    tiles = recurse(tiles, neighborTileNumber, maximumRecursions, ...
                    recursionNumber);

  % if neighbor tile is already in a different cluster,
  elseif tiles(neighborTileNumber).clusterNumber ~= ...
         tiles(tileNumber).clusterNumber,
    
    % merge current cluster into the other cluster
    for mergeTileNumber = 1 : numberOfTiles,
      if tiles(mergeTileNumber).clusterNumber == ...
         tiles(tileNumber).clusterNumber,
        tiles(mergeTileNumber).clusterNumber = ...
            tiles(neighborTileNumber).clusterNumber;
      end
    end

  % if neighbor tile is already in this cluster,
  else

  % continue
  end

% end loop over neighboring tiles
end

% return to calling function
return;
