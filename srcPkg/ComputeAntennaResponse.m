function [Fp, Fc] = ComputeAntennaResponse(phi,theta,psi,detector)
% ComputeAntennaResponse - Compute the antenna response factors for a
% specfic detector and set of sky positions and polarizations.  
%
%   [Fp, Fc] = ComputeAntennaResponse(phi,theta,psi,detector)
%
%  phi      Vector.  Azimuthal angle of the sky position of the source in
%           Earth-centered coordinates.  The prime meridian at Greenwich is
%           phi=0.
%  theta    Vector.  Polar angle of the sky position of the source in
%           Earth-centered coordinates.  The north pole is at theta=0; the
%           south pole is at theta=pi.
%  psi      Vector.  The polarization angle of the gravitational wave,
%           defined as the angle counterclockwise about the direction of
%           PROPAGATION from a line of constant theta pointing to
%           decreasing phi to the positive x axis of the source coordinates
%           (the source "+" polarization direction).
%  detector String or 3x3 numerical array.  If a string, it specifies the 
%           detector site on to which to project the signal, and must be 
%           one of the values recognized by the function LoadDetectorData.
%           If an array, then it specifies the detector response matrix 
%           d^{ab} to be used.
%
%   Fp      Column vector.  The "plus" antenna response factors for the  
%           specified detector.
%   Fc      Column vector.  The "cross" antenna response factors for the  
%           specified detector.
%
% The vectors phi, theta, psi must have the same length, except that one or
% more may be scalar, in which case the scalars are expanded to vectors of
% the same size.
%
% initial write: Patrick J. Sutton 2004.04.26
%
% $Id: ComputeAntennaResponse.m 1941 2007-11-14 05:00:57Z jrollins $

%----- Make sure input arguments are column vectors.
if (size(phi,2)>1)
    phi = phi(:);
end
if (size(theta,2)>1)
    theta = theta(:);
end
if (size(psi,2)>1)
    psi = psi(:);
end

%----- Verify that input angle vectors are the same size, or are scalars.
maxLength = max([length(theta), length(phi), length(psi)]);
if (length(phi)~=maxLength)
    % ---- Make sure it's a scalar, then expand to vector.
    if (isscalar(phi))
        phi = phi*ones(maxLength,1);
    else
        error(['Sky position, polarization angles phi, theta, psi must ' ...
            'have the same length.'])
    end
end
if (length(theta)~=maxLength)
    % ---- Make sure it's a scalar, then expand to vector.
    if (isscalar(theta))
        theta = theta*ones(maxLength,1);
    else
        error(['Sky position, polarization angles phi, theta, psi must ' ...
            'have the same length.'])
    end
end
if (length(psi)~=maxLength)
    % ---- Make sure it's a scalar, then expand to vector.
    if (isscalar(psi))
        psi = psi*ones(maxLength,1);
    else
        error(['Sky position, polarization angles phi, theta, psi must ' ...
            'have the same length.'])
    end
end

%----- Get the detector response matrix "d^{ab}_i".
%      If input argument is a detector/site name, retrieve this data by
%      calling LoadDetectorData.  If the input argument is a 3x3 numerical 
%      array, use that for the detector response matrix. 
if (ischar(detector))
   %----- Load data on detector.
   DetData = LoadDetectorData(detector);
   d = DetData.d;
elseif (isnumeric(detector) && isequal(size(detector),[3 3]))
   d = detector;
else
   error(['Detector not recognized. 4th argument should be a ' ... 
       'detector/site name or a 3x3 array.']);  
end
%----- Convert to vector.
d = d(:);

%----- Compute polarization tensors (functions of the sky position and
%      polarization). 
m1 = [ sin(phi).*cos(psi)-cos(phi).*cos(theta).*sin(psi) ];
m2 = [-cos(phi).*cos(psi)-sin(phi).*cos(theta).*sin(psi) ];
m3 = [ sin(theta).*sin(psi) ];
n1 = [-sin(phi).*sin(psi)-cos(phi).*cos(theta).*cos(psi) ];
n2 = [ cos(phi).*sin(psi)-sin(phi).*cos(theta).*cos(psi) ];
n3 = [ sin(theta).*cos(psi) ];
mm = [ m1.*m1, m1.*m2, m1.*m3, m2.*m1, m2.*m2, m2.*m3, m3.*m1, m3.*m2, m3.*m3 ];
mn = [ m1.*n1, m1.*n2, m1.*n3, m2.*n1, m2.*n2, m2.*n3, m3.*n1, m3.*n2, m3.*n3 ];
nm = [ n1.*m1, n1.*m2, n1.*m3, n2.*m1, n2.*m2, n2.*m3, n3.*m1, n3.*m2, n3.*m3 ];
nn = [ n1.*n1, n1.*n2, n1.*n3, n2.*n1, n2.*n2, n2.*n3, n3.*n1, n3.*n2, n3.*n3 ];
e_plus = mm - nn;
e_cross = mn + nm;

%----- Compute waveform projected onto antenna pattern.
Fp = e_plus*d;
Fc = e_cross*d;

% % ---- This direct matrix form is clearer, but only works for a single
% %      sky position and polarization angle at a time.
% m = [  sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi) ; ...
%       -cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi) ; ...
%        sin(theta)*sin(psi) ];
% n = [ -sin(phi)*sin(psi)-cos(phi)*cos(theta)*cos(psi) ; ...
%        cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi) ; ...
%        sin(theta)*cos(psi) ];
% e_plus = m*m' - n*n';
% e_cross = m*n' + n*m';
% Fp = trace(e_plus*d);
% Fc = trace(e_cross*d);

%----- Done
return
