function gmst = gpstogmst(gpst)
% GPSTOGMST Convert from GPST to GMST
%
% GPSTOGMST converts the given Global Positioning System Time (GPST) to
% Greenwich Mean Sidereal Time (GMST).
%
% usage:  gmst = gpstogmst(gpst)
%
%   gpst    Global Positioning System Time (GPST) [seconds]
%
%   gmst    Greenwich Mean Sidreal Time (GMST) [rad]
%
% Algorithm taken from the International Astronomical Union (IAU)
% Standards Of Fundamental Astronomy (SOFA) Libraries from 2008
% March 01.  Uses the GMST82 algorithm based on IAU 1982 model to
% convert from UT1 to GMST.
% 
% http://www.iau-sofa.rl.ac.uk/2008_0301/sofa/gmst82.html
%
% UTC is used as an estimate of UT1.  GPST is converted to UTC by
% removing leap seconds announced by the International Earth Rotation
% and Reference Systems Service (IERS).
% 
% http://hpiers.obspm.fr/eop-pc/
%
% This list of leap seconds must be updated manually when a new leap
% second is announced.  The most recent update is for the leap second
% at 2008-12-31 12:59:60 UTC.

% Shourov K. Chatterji <shourov@ligo.caltech.edu>

% %Id%

% list of leap seconds since 1980-01-01
leapData = [046828800 +1;...
	    078364801 +1;...
	    109900802 +1;...
	    173059203 +1;...
	    252028804 +1;...
	    315187205 +1;...
	    346723206 +1;...
	    393984007 +1;...
	    425520008 +1;...
	    457056009 +1;...
	    504489610 +1;...
	    551750411 +1;...
	    599184012 +1;...
	    820108813 +1;...
	    914803214 +1;];

% limit of validity of leap second data
leapValid = 930441614;

% number of leap seconds
numberOfLeaps = size(leapData, 1);

% reference GPS time (2000-01-01 12:00:00 UTC) for GMST conversion
reference = 630763213;

% issue warning if outside of valid leap second data
if any(gpst > leapValid),
  warning('leap seconds unknown');
end

% determine leap second adjustment
leapSeconds = zeros(size(gpst));
for leapNumber = 1 : numberOfLeaps,
  if leapData(leapNumber, 1) >= reference,
    leapSeconds = leapSeconds + ...
        leapData(leapNumber, 2) * (leapData(leapNumber, 1) < gpst);
  else
    leapSeconds = leapSeconds - ...
        leapData(leapNumber, 2) * (leapData(leapNumber, 1) >= gpst);
  end
end

% universal time since reference epoch
seconds = gpst - reference - leapSeconds;
days = seconds / 86400;
centuries = days / 36525;

% fractional day since zero hour universal time
fraction = mod(days, 1) - 0.5;

% greenwich mean sidereal time in seconds
gmst = +2.411054841000e+04 + ...
       +8.640184812866e+06 * centuries + ...
       +9.310400000000e-02 * centuries.^2 + ...
       -6.200000000000e-06 * centuries.^3 + ...
       86400 * fraction;

% restrict to [0, 86400) seconds
gmst = mod(gmst, 86400);

% convert to radians
gmst = 2 * pi * gmst / 86400;

% return
return
