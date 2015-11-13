function distances = wdistance(significants, distanceMetric, ...
                               durationInflation, bandwidthInflation)
% WDISTANCE Compute distances between significant Q transform tiles
%
% WDISTANCE computes the distance between significant Q transform
% tiles produced by WTHREHSOLD or WSELECT.
%
% distances = wdistance(significants, distanceMetric, ...
%                       durationInflation, bandwidthInflation);
%
%   distances            cell array of significant tiles distances
%
%   significants         cell array of significant tiles properties
%   distanceMetric       choice of metric for computing distance
%   durationInflation    multiplicative scale factor for duration
%   bandwidthInflation   multiplicative scale factor for bandwidth
%
% WDISTANCE expects a cell array of Q transform event structures with
% one cell per channel.  The event structures must contain at least
% the following fields, which describe the properties of statistically
% significant tiles used to compute distance.
% 
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   q                    quality factor of tile []
%   normalizedEnergy     normalized energy of tile []
%
% WDISTANCE returns a cell array of Q transform distance structures
% with one cell per cahnnel.  In addition to a structure identifier,
% the distance structures contains the following single field.
%
%   distance             pairwise distance between tiles
%
% Distances are returned in the same format as the PDIST function.
% In particular, distances are reported as row vectors of length
% N * (N - 1) / 2, where N is the number of significant tiles for
% a given channel.  This row vector is arranged in the order of
% (1,2), (1,3), ..., (1,N), (2,3), ..., (2,N), ..., (N-1, N).  Use
% the SQUAREFORM function to convert distances into a matrix format.
%
% The following choices of distance metric are provided
%
%   'pointMismatch'
%
%     The second order expansions of the mismatch metric used for
%     tiling the signal space, evaluated at the center point between
%     two tiles.
%
%   'integratedMismatch'
%
%     The integrated second order expansion of the mismatch metric
%     between two tiles.  The center point between two tiles is used
%     to determine the metric coefficients.  Only the diagonal metric
%     terms are integrated, and these are added in quadrature.
%
%   'logMismatch'
%
%     The exact mismatch between two tiles.  The result is returned
%     as the negative natural logarithm of the overlap between two
%     tiles.
%
%   'euclidean'
% 
%     The dimensionless euclidean distance between tiles in the
%     time-frequency plane after normalizing by the mean duration and
%     bandwidth of each tile pair.
%
%   'modifiedEuclidean'
% 
%     A modified calculation of euclidean distance based on emperical
%     studies that gives more weight to differences in frequency than
%     in time.
%
% The optional durationInflation and bandwidthInflation arguments are
% multiplicative scale factors that are applied to the duration and
% bandwidth of significant tiles prior to determining their distance.
% If not specified, these parameters both default to unity such that
% the resulting tiles have unity time-frequency area.  They are only
% used in the calculation of the euclidean metric distance.
%
% See also WTHRESHOLD, WSELECT, WCLUSTER, SQUAREFORM, and PDIST.

% Rubab Khan
% rmk2109@columbia.edu

% Shourov Chatterji
% shourov@ligo.caltech.edu

% $Id: wdistance.m 986 2008-08-14 21:04:56Z lstein $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 4, nargin));

% apply default arguments
if (nargin < 2) || isempty(distanceMetric),
  distanceMetric = 'integratedMismatch';
end
if (nargin < 3) || isempty(durationInflation),
  durationInflation = 1.0;
end
if (nargin < 4) || isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
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
%                       initialize distances structures                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of significant tile distances
distances = cell(size(significants));

% begin loop over channels
for channelNumber = 1 : numberOfChannels

  % insert structure identification string
  distances{channelNumber}.id = 'Discrete Q-transform distance structure';
  
% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      create pairwise list of tiles                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of significant tiles
  numberOfSignificants = length(significants{channelNumber}.time);

  % number of unique significant tile pairs
  numberOfPairs = numberOfSignificants * (numberOfSignificants - 1) / 2;

  % build list of pairwise indices
  [pairIndices1, pairIndices2] = meshgrid(1 : numberOfSignificants);
  pairIndices1 = tril(pairIndices1, -1);
  pairIndices2 = tril(pairIndices2, -1);
  pairIndices1 = pairIndices1(:);
  pairIndices2 = pairIndices2(:);
  pairIndices1 = pairIndices1(find(pairIndices1));
  pairIndices2 = pairIndices2(find(pairIndices2));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                        determine tile properties                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % extract significant tile properties
  time1 = significants{channelNumber}.time(pairIndices1);
  time2 = significants{channelNumber}.time(pairIndices2);
  frequency1 = significants{channelNumber}.frequency(pairIndices1);
  frequency2 = significants{channelNumber}.frequency(pairIndices2);
  q1 = significants{channelNumber}.q(pairIndices1);
  q2 = significants{channelNumber}.q(pairIndices2);
  normalizedEnergy1 = significants{channelNumber}.normalizedEnergy(pairIndices1);
  normalizedEnergy2 = significants{channelNumber}.normalizedEnergy(pairIndices2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         determine tile distances                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % switch on distance metric
  switch lower(distanceMetric),
  
    % handle case of point mismatch metric distance
    case 'pointmismatch',
      
      % mean properties of pairs of tiles
      meanFrequency = sqrt(frequency1 .* frequency2);
      meanQ = sqrt(q1 .* q2);

      % parameter distances between pairs of tiles
      timeDistance = time2 - time1;
      frequencyDistance = frequency2 - frequency1;
      qDistance = q2 - q1;

      % mismatch metric distances between pairs of tiles
      distances{channelNumber}.distance = ...
          timeDistance.^2 * 4 * pi^2 .* meanFrequency.^2 ./ meanQ.^2  + ...
          frequencyDistance.^2 .* (2 + meanQ.^2) ./ (4 * meanFrequency.^2)  + ...
          qDistance.^2 ./ (2 * meanQ.^2) - ...
          frequencyDistance .* qDistance ./ (meanFrequency .* meanQ);
  
    % handle case of integrated mismatch metric distance
    case 'integratedmismatch',
      
      % mean properties of pairs of tiles
      meanFrequency = sqrt(frequency1 .* frequency2);
      meanQ = sqrt(q1 .* q2);

      % parameter distances between pairs of tiles
      timeDistance = 2 * pi * meanFrequency .* (time2 - time1) ./ meanQ;
      frequencyDistance = sqrt(2 + meanQ.^2) .* log(frequency2 ./ frequency1) / 2;
      qDistance = log(q2 ./ q1) / sqrt(2);

      % mismatch metric distances between pairs of tiles
      distances{channelNumber}.distance = ...
          sqrt(timeDistance.^2 + frequencyDistance.^2 + qDistance.^2);
  
    % handle case of log mismatch distance
    case 'logmismatch',
    
      % report error
      error('logMismatch metric not yet implemented');
  
    % handle case of euclidean distance
    case 'euclidean',
    
      % tile dimensions
      bandwidth1 = 2 * sqrt(pi) * frequency1 ./ q1;
      bandwidth2 = 2 * sqrt(pi) * frequency2 ./ q2;
      duration1 = 1 ./ bandwidth1;
      duration2 = 1 ./ bandwidth2;

      % apply tile inflation factors
      duration1 = duration1 * durationInflation;
      duration2 = duration2 * durationInflation;
      bandwidth1 = bandwidth1 * bandwidthInflation;
      bandwidth2 = bandwidth2 * bandwidthInflation;

      % time and frequency distance scales
      timeScale = (duration1 .* normalizedEnergy1 + ...
                   duration2 .* normalizedEnergy2) / ...
                  (normalizedEnergy1 + normalizedEnergy2);
      frequencyScale = (bandwidth1 .* normalizedEnergy1 + ...
                        bandwidth2 .* normalizedEnergy2) / ...
                       (normalizedEnergy1 + normalizedEnergy2);

      % normalized time and frequency distance between tiles
      timeDistance = abs(time2 - time1) / timeScale;
      frequencyDistance = abs(frequency2 - frequency1) / frequencyScale;

      % normalized euclidean distance between tiles
      distances{channelNumber}.distance = ...
          sqrt(timeDistance.^2 + frequencyDistance.^2);
  
    % handle case of modified euclidean distance
    case 'modifiedeuclidean',
    
      % tile dimensions
      bandwidth1 = 2 * sqrt(pi) * frequency1 ./ q1;
      bandwidth2 = 2 * sqrt(pi) * frequency2 ./ q2;
      duration1 = 1 ./ bandwidth1;
      duration2 = 1 ./ bandwidth2;

      % apply tile inflation factors
      duration1 = duration1 * durationInflation;
      duration2 = duration2 * durationInflation;
      bandwidth1 = bandwidth1 * bandwidthInflation;
      bandwidth2 = bandwidth2 * bandwidthInflation;

      % time and frequency distance scales
      timeScale = (duration1 .* normalizedEnergy1 + ...
                   duration2 .* normalizedEnergy2) / ...
                  (normalizedEnergy1 + normalizedEnergy2);
      frequencyScale = (bandwidth1 .* normalizedEnergy1 + ...
                        bandwidth2 .* normalizedEnergy2) / ...
                       (normalizedEnergy1 + normalizedEnergy2);

      % normalized time and frequency distance between tiles
      timeDistance = abs(time2 - time1) / timeScale;
      frequencyDistance = abs(frequency2 - frequency1) / frequencyScale;

      % modified normalized euclidean distance between tiles
      distances{channelNumber}.distance = ...
          sqrt(timeDistance.^2 + 30 * frequencyDistance.^2);

    % handle unknown distance metric
    otherwise,

      % report error
      error(['unknown distance metric "' distanceMetric '"']);
      
  % end switch on distance metric
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   return statistically significant events                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
