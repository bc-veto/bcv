function coordinates = wconvertskycoordinates(coordinates, time, system0, system1)
% WCONVERTSKYCOORDINATES convert between sky coordinate systems
%
% WCONVERTSKYCOORDINATES for given input coordinates in one
% coordinate system, convert to a different coordinate system
% 
% usage:
%   coordinates = wconvertskycoordinates(coordinates, time, system0, system1)
%
% coordinates must be a 2-column matrix of coordinates in the given system.
% time maybe be a scalar, or a vector of the same length as coordinates.
% Available systems for conversion are:
%
%   'geocentric'  [theta, phi]  (radians)
%   'equatorial'  [dec, ra]     (degrees)
%
% See also WTILESKY, WRESPONSE, and WSKYMAP.

% Jameson Rollins <jrollins@phys.columbia.edu>
% Patrick J. Sutton <Patrick.Sutton@astro.cf.ac.uk>

% $Id$

possibleSystems = {'geocentric','equatorial'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(4, 4, nargin));

if size(coordinates,2) ~= 2,
    error('coordinates must have two columns');
end
if ~any(strcmpi(system0, possibleSystems)),
  error('unrecognized input coordinate system');
end

if ~any(strcmpi(system1, possibleSystems)),
  error('unrecognized output coordinate system');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              convert coordinates                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = [system0,'2',system1];
coordinates = feval(fname,coordinates,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return to calling function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                subfunctions                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             equatorial2geocentric                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_phi = equatorial2geocentric(dec_ra,gps)
% EQUATORIAL2GEOCENTRIC Convert from equatorial to geocentric coordinates.
%
% usage:
%   theta_phi = equatorial2geocentric(dec_ra,gps)
%
% dec_ra     2-column matrix of declinations and right accensions (degrees)
% gps        Vector or scalar of GPS times (seconds).  If a scalar, the same
%            GPS time is used for converting each sky position.
%
% theta_phi  2-column vector of polar and azimuthal angles (radians)
%
% Geocentric coordinates have the x-axis (phi=0) pointing from the
% Earth's center to the intersection of the prime meridian of
% Greenwich with the equator, the z-axis (theta=0) pointing from the
% Earth's center to the North pole, and the y-xais chosen to form a
% right-handed coordinate system.

% Convert everything to column vectors.
dec = dec_ra(:,1);
ra = dec_ra(:,2);
gps = gps(:);

% check coordinates
if any((dec < -90) | (dec > 90)),
  error('declination outside of range [-90, 90]');
end
if any((ra < 0) | (ra >= 360)),
  error('right accension outside of range [0, 360)');
end

% Compute sidereal time (sec) at each event time.
gmst = gps2gmst(gps);
% Convert to degrees
gmstDeg = gmst / 86400 * 360;

% Convert to vector if needed.
if length(gmstDeg)==1
    gmstDeg = gmstDeg * ones(size(dec));
end

% Compute geocentric coordinates, in degrees.
thetaDeg  = 90 - dec;
phiDeg = ra - gmstDeg;

% Convert to radians.
theta = thetaDeg * pi / 180;
phi = phiDeg * pi / 180;
phi = mod(phi, 2 * pi);

theta_phi = [theta phi];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             equatorial2equatorial                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dec_ra = equatorial2equatorial(dec_ra,gps)
% EQUATORIAL2EQUATORIAL trivial convertion to same

% Convert everything to column vectors.
dec = dec_ra(:,1);
ra = dec_ra(:,2);

% check coordinates
if any((dec < -90) | (dec > 90)),
  error('declination outside of range [-90, 90]');
end
if any((ra < 0) | (ra >= 360)),
  error('right accension outside of range [0, 360)');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             geocentric2equatorial                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dec_ra = geocentric2equatorial(theta_phi,gps)
% GEOCENTRIC2EQUATORIAL Convert from geocentric to equatorial coordinates.
%
% usage:
%   dec_ra = equatorial2geocentric(theta_phi,gps)
%
% theta_phi  2-column vector of polar and azimuthal angles (radians)
% gps        Vector or scalar of GPS times (seconds).  If a scalar, the same
%            GPS time is used for converting each sky position.
%
% dec_ra     2-column matrix of declinations and right accensions (degrees)
%
% Geocentric coordinates have the x-axis (phi=0) pointing from the
% Earth's center to the intersection of the prime meridian of
% Greenwich with the equator, the z-axis (theta=0) pointing from the
% Earth's center to the North pole, and the y-xais chosen to form a
% right-handed coordinate system.

% Convert everything to column vectors.
theta = theta_phi(:,1);
phi = theta_phi(:,2);
gps = gps(:);

% check coordinates
if any((theta < 0) | (theta > pi)),
  error('theta outside of range [0, pi]');
end
if any((phi < 0) | (phi >= 2 * pi)),
  error('phi outside of range [0 , 2 pi)');
end

% Compute sidereal time (sec) at each event time.
gmst = gps2gmst(gps);
% Convert to degrees
gmstDeg = gmst / 86400 * 360;

% Convert to vector if needed.
if length(gmstDeg)==1
  gmstDeg = gmstDeg * ones(size(theta));
end

% Convert to degrees
thetaDeg = theta * 180 / pi;
phiDeg = phi * 180 / pi;

% Compute equatorial coordinates, in degrees.
dec = 90 - thetaDeg;
ra = phiDeg + gmstDeg;
ra =  mod(ra,360);

dec_ra = [dec ra];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             geocentric2geocentric                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_phi = geocentric2geocentric(theta_phi,gps)
% GEOCENTRIC2GEOCENTRIC trivial convertion to same

% Convert everything to column vectors.
theta = theta_phi(:,1);
phi = theta_phi(:,2);

% check coordinates
if any((theta < 0) | (theta > pi)),
  error('theta outside of range [0, pi]');
end
if any((phi < 0) | (phi >= 2 * pi)),
  error('phi outside of range [0, 2 pi)');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    gps2gmst                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gmst, gmstHms] = gps2gmst(gps)
% GPS2GMST Convert from GPS time to Greenwich Mean Sidereal Time.
%
% usage:
%   [gmst, gmstHms] = gps2gmst(gps)
%
% gps       Input GPS time.  May be a scalar or a vector.  If a vector,
%           then the output variables will be arrays with one row per gps
%           value.
%
% gmst      Greenwich Mean Sidereal Time (s)
% gmstHms   Greenwich Mean Sidereal Time (h:m:s)

% Note: If gps is a vector then gmstH, etc will have same dimensions.
% Make sure we covert to column vectors before assembling gmstHms array.
if (min(size(gps))==0)
    gmst = [];
    gmstHms = [];
    return;
elseif (min(size(gps))==1)
    if (size(gps,2)>1)
        gps = gps';
    end
else
    error('Input time must be a scalar or vector.');
end

% GPS time of J2000 epoch (approx 2000 Jan 12, 12:00 UTC)
gps0 = 630763213;

% (float) Days since J2000
D = (gps-gps0)/86400;

% (int+1/2) Days between J2000 and last 0h at Greenwich (always integer + 1/2)
d_u = floor(D)+0.5;
d_u(logical((D-floor(D))<0.5)) = d_u(logical((D-floor(D))<0.5)) - 1;

% (float) Fraction of day since last 0h at Greenwich
df_u = D-d_u;

% GMST (s) at Greenwich at last 0h:
T_u = d_u/36525;
gmst0h = 24110.54841 + 8640184.812866*T_u + 0.093104*T_u.^2 - 6.2e-6*T_u.^3;

% Current GMST (s)
gmst = gmst0h + 1.00273790935*86400*df_u;

% Remove any integer days
if (gmst>=86400)
    gmst = gmst - floor(gmst/86400)*86400;
end

% In hours:min:sec
gmstH = floor(gmst/3600);
gmstM = floor((gmst-gmstH*3600)/60);
gmstS = (gmst-gmstH*3600-gmstM*60);
gmstHms = [gmstH gmstM gmstS];

return
