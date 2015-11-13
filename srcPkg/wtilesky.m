function [coordinates, solidAngles, probabilities] = ...
        wtilesky(angularResolution, input, coordinateSystem)
% WTILESKY Tile the sky with a specified angular resolution
%
% WTILESKY returns a set of spherical coordinates that covers the sky, or a
% specified region of the sky, with close to the minimum number of points
% required to ensure an angular resolution consistent with the frequency
% range of the search.
%
% usage:
%   [coordinates, solidAngles, probabilities] = ...
%             wtilesky(angularResolution, coordinates, coordinateSystem);
%
% all angular measurements should be given in radians.
%
%   angularResolution  angular resolution to use in tiling
%   coordinates        specification of sky coordinates
%   coordinateSystem   coordinate system of sky coordinates
%
% If angularResolution is a structure, it is assumed to be a Q-transform tiling
% structure from wtile, and the resolution is determined from the highest
% frequency in the tile.
%
% The 'coordinates' can be given in a number of different ways:
%
%   'all'              tile the entire sky (default)
%   N                  tile entire sky with N points (ignore resolution)
%   [a b]              2-column matrix of sky coordinates
%   {[aMin aMax] [bMin bMax]}
%                      2-element cell of 2-element vectors that specify ranges
%   'filename'         filename containing a 2-column matrix of sky coordinates
%
% The 'coordinateSystem' choices, with their labels, are:
%
%   'geocentric'       [theta phi] (default)
%   'equatorial'       [dec ra]
%
% The outputs are:
%
%   coordinates        matrix of sky coordinates (radians)
%   solidAngles        vector of solid angles (steradians)
%   probabilities      vector of a priori uniform probabilities
%
% The coordinate theta is a geocentric colatitude running from [0, pi],
% with 0 at the North pole to pi at the South pole.  The coordinate
% phi is the geocentric longitude running from [0, 2 pi), with 0 on the
% prime meridian.
%
% WTILESKY also returns the approximate differential solid angle associated with
% each sky coordinates, and its corresponding a priori uniform probability density.
%
% The angular resolution of the sky tiling is determined from the specified q
% transform tiling structure in order to ensure a worst case phase mismatch of
% pi/4 for signals at the maximum frequency of the analysis.
%
% WTILESKY uses the HEALPix algorithm (http://healpix.jpl.nasa.gov/).
%
% See also WCONVERTSKYCOORDINATES, WTILE, WRESPONSE, and WSKYMAP.

% Shourov Chatterji <shourov@ligo.caltech.edu>
% Albert Lazzarini <lazz@ligo.caltech.edu>
% Antony Searle <antony.searle@anu.edu.au>
% Leo Stein <lstein@ligo.caltech.edu>
% Patrick Sutton <psutton@ligo.caltech.edu>
% Massimo Tinto <massimo.tinto@jpl.nasa.gov>
% Yoichi Aso <aso@astro.columbia.edu>
% Jameson Rollins <jrollins@phys.columbia.edu>

% $Id: wtilesky.m 909 2008-07-02 10:19:07Z shourov $

possibleSystems = {'geocentric','equatorial'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 4, nargin));

% apply default arguments
if (nargin < 2) || isempty(input),
    input = 'all';
end
if (nargin < 3) || isempty(coordinateSystem),
    coordinateSystem = 'geocentric';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             check coordinate system                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~any(strcmpi(coordinateSystem,possibleSystems)),
  error(['unknown coordinate system "',coordinateSystem,'"'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               check coordinates                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordinates = [];
range = [];
points = [];

% check coordinates specification
switch class(input)
  case 'char'
    switch lower(input)
      case 'all',
        % leave variables empty for all sky
      otherwise
        % assume file name and attempt to load positions from file
        coordinates = load(input);
    end
  case 'double'
    if length(input) == 1,
      % input is number of points in tile
      points = input;
    else
      % vector of sky coordinatess (theta, phi)
      coordinates = input;
    end
  case 'cell'
    % 2-element cell array with range of first and second coordinate
    range = input;
  otherwise
    error('unknown coordinates specification')
end

% check coordinate or range values
if ~isempty(coordinates),

  if size(coordinates,2) ~= 2,
    error('coordinate must be scalar or two column vector')
  end

  % check that coordinates fall in appropriate range for coordinate system
  feval(['checkRange_',coordinateSystem],coordinates(:,1),coordinates(:,2));

elseif ~isempty(range),

  if length(range{1}) ~= 2 || length(range{2}) ~= 2,
    error('coordinate ranges must be two element vectors [min max]');
  end

  % check that ranges fall in appropriate range for coordiante system
  feval(['checkRange_',coordinateSystem],range{1},range{2});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            set angular resolution                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if the points variable is not empty
if ~isempty(points),

  % set the angular resolution to be whole sky (4pi) divided by number of points
  angularResolution = sqrt(4*pi/points);

% otherwise
else

  % if it's a structure, determine resolution from max Q tile frequency
  if isstruct(angularResolution),
    tiling = angularResolution;
    % validate tiling structure
    if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
        error('input argument is not a discrete Q transform tiling structure');
    end

    % worst case phase shift due to wrong sky coordinates
    phaseResolution = pi / 4;

    % required timing resolution at maximum frequency of search
    timingResolution = phaseResolution / (2 * pi * tiling.lowPassCutoff);

    % maximum possible time delay between earth based detectors
    earthRadius = 6378135;
    speedOfLight = 299792458;
    maximumTimeDelay = 2 * earthRadius / speedOfLight;

    % required angular resolution
    angularResolution = pi * timingResolution / maximumTimeDelay;

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            healpix tiling of sky                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solidAngles = [];
probabilities = [];

% tile sky with geocentric coordinates
if isempty(coordinates),

  % The tiling is specified by a parameter N, which is required to be a
  % power of 2.  The number of pixels on the sky is 12*N^2.  The solid
  % angle of each pixel is pi/(3*N^2).  At the moment, the resolution is
  % then assumed (with some inaccuracy) smallest power of 2 bigger than
  % the square-root of the solid angle of the pixels.  This should be
  % done more accurately based on the maximum angular separation between
  % pixels in the tiling.
  %
  % General strategy for the coding
  % In HEALPix tiling the sky is divided into four parts, north polar
  % cap, north equatorial band, south equatorial band, south polar
  % cap. For the definition of those regions, please refer to the
  % reference above.  Since the formula for the index-to-coordinates
  % conversion differs from region to region, we divide the indices
  % vector p into 4 parts according to the values.

  % determine total number of points from angularResolution
  Npix = 4*pi/(angularResolution^2);

  % determine resolution parameter from angular resolution
  % N = N_side from paper
  N = sqrt(Npix/12);
  N = 2^(ceil(log2(N)));
  Npix = 12*N^2;

  % Prepare a sky position index vector
  p=[0:Npix-1];
  p=p(:);

  % Divide the indices p into 4 parts according to the values.
  % North Polar cap
  indices1=find(p<2*N*(N-1));
  p1=p(indices1);  % Extract p in this region
  ph=(p1+1)/2;
  i=floor(sqrt(ph-sqrt(floor(ph))))+1;
  j=p1+1-2*i.*(i-1);
  z1 = 1-(i.^2)/(3*N^2);
  phi1 = pi./(2*i).*(j-1/2);

  % North Equatorial belt
  indices2=find((p>=2*N*(N-1)).*(p<Npix/2));
  p2=p(indices2);
  ph=p2-2*N*(N-1);
  i=floor(ph/(4*N))+N;
  j=mod(ph, 4*N)+1;
  s=mod(i-N+1, 2);
  z2 = 4/3 - 2*i/(3*N);
  phi2 = pi./(2*N).*(j-s/2);

  % South Equatorial belt
  z3 = -z2;
  phi3 = -phi2;

  % South Polar cap
  z4 = -z1;
  phi4 = -phi1;

  % Combine the results from each region
  z = [z1;z2;z3;z4];
  theta = acos(z);

  phi = [phi1;phi2;phi3;phi4];

  % force phi in range [0, 2 pi)
  phi = mod(phi, 2*pi);

  % Solid angle (it is constant).
  solidAngles = pi/(3*N^2)*ones(size(theta));
  % A priori probability
  probabilities = 1 / (4*pi) * solidAngles;

  coordinates = [theta(:) phi(:)];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                           convert coordinates                              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % convert to desired coordinate system if not geocentric
  if ~strcmpi(coordinateSystem,'geocentric'),

      % convert to geocentric coordinates, with dummy time of 0,
      % since we're not actually doing a real conversion.
      coordinates = wconvertskycoordinates(coordinates, 0, ...
                                           'geocentric', coordinateSystem);

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            filter coordinate range                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(range),

  % output coordinates only in specified desired range
  coordinates = coordinates(find(...
      coordinates(:,1) >= range{1}(1) & coordinates(:,1) <= range{1}(2) & ...
      coordinates(:,2) >= range{2}(1) & coordinates(:,2) <= range{2}(2)),:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        solid angles and probabilities                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate coordinate solid angles if haven't already
if isempty(solidAngles),

  theta = coordinates(:,1);
  solidAngles = pi/(3*angularResolution^2)*ones(size(theta));

end

% calculate a priori probabilities if haven't already
if isempty(probabilities),

  probabilities = 1 / (4*pi) * solidAngles;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                end main block                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                subfunctions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          geocentric coordinate check                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkRange_geocentric(theta,phi)
% first coordinate: theta
if any((theta < 0) | (theta > pi)),
   error('theta outside of range [0, pi]');
end

% second coordinate: phi
if any((phi < 0) | (phi >= 2 * pi)),
  error('phi outside of range [0, 2 pi)');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         equatorial coordinate check                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkRange_equatorial(dec,ra)
% first coordinate: declination
if any((dec < -90) | (dec > 90)),
  error('declination outside of range [-90 90]');
end

% first coordinate: right accension
if any((ra < 0) | (ra >= 360)),
  error('right accension outside of range [0 360)');
end

return
