function baselines = wbaseline(detectors)
% WBASELINE Light travel time between detectors
%
% WBASELINE computes the matrix of pairwise distances between specified
% gravitational wave detectors in terms of the light travel time.
%
%  usage:
%
%    baselines = wbaseline(detectors);
%
%    detectors      cell array of site, detector, or channel names
%
%    baselines      matrix of pairwise detector distances [seconds]
%
% The list of detectors should be specified as a cell array of site identifiers
% or channel names, which are passed to the WDETECTORS function to retrieve the
% corresponding detector location and response data.  Site identifiers are
% determined by taking the first character of each channel name.  The following
% sites are recognized by the WDETECTORS function.
%
%   G  GEO
%   H  LIGO Hanford
%   L  LIGO Livingston
%   V  Virgo
%   0  Fiducial reference detecor
%
% The fiducial reference detector is provided for validation purposes.  It is
% positioned in the direction of the North pole, at a distance equal to one
% speed of light travel time from the Earth center, and has arms along the
% positive X and Y axis in the Earth fixed geocentric coordinate system.
%
% The elements of the resulting baseline matrix are the speed of light travel
% times between pairs of detector sites in units of seconds.  The matrix is
% indexed in the same order as the input cell array of detector names.
%
% See also WDETECTORS, WRESPONSE, WTILESKY, and WSKYMAP.
%
% Authors:
% Shourov Chatterji <shourov@ligo.caltech.edu>
% Albert Lazzarini <lazz@ligo.caltech.edu>
% Antony Searle <antony.searle@anu.edu.au>
% Leo Stein <lstein@ligo.caltech.edu>
% Patrick Sutton <psutton@ligo.caltech.edu>
% Massimo Tinto <massimo.tinto@jpl.nasa.gov>

% $Id: wbaseline.m 678 2008-04-18 19:21:39Z acsearle $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 1, nargin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% force cell array of detector names
detectors = wmat2cell(detectors);

% force one dimensional cell array
detectors = detectors(:);

% number of detectors
numberOfDetectors = length(detectors);

% force uppercase detector names
detectors = upper(detectors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              detector geometry                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% retrieve detector location data
locations = wdetectors(detectors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          compute baseline distances                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize results
baselines = zeros(numberOfDetectors, numberOfDetectors);

% construct unique distances
for detectorNumber1 = 1 : numberOfDetectors,
  for detectorNumber2 = detectorNumber1 + 1 : numberOfDetectors,
    baselines(detectorNumber1, detectorNumber2) = ...
      norm(locations{detectorNumber1} - locations{detectorNumber2});
  end
end

% create symmetric baseline matrix
baselines = baselines + baselines.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return to calling function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
