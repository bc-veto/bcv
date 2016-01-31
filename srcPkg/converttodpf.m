function [Fp, Fc, psi] = converttodpf(Fp,Fc)
%
% convertToDominantPolarizationFrame - Convert network antenna response 
% vectors to dominant polarization frame.
%
%   [FpDP FcDP psiDP] = converttodpf(Fp,Fc)
%
%   Fp      Array of "plus" antenna responses.  Each row holds the antenna
%           responses for one sky position.
%   Fc      Array of "cross" antenna responses.  Each row holds the antenna
%           responses for one sky position.
%
%   FpDP    Array of "plus" antenna responses in the dominant polarization
%           frame. 
%   FcDP    Array of "cross" antenna responses in the dominant polarization
%           frame. 
%   psiDP   Angle between the input polarization frame and the dominant
%           polarization frame, computed for each sky position.
%
% The dominant polarization frame is defined as the frame in which the
% network antenna repsonse vectors Fp "plus" and Fc "cross" are orthogonal
% and |Fp| > |Fc|.
%
% initial write: Patrick J. Sutton 2006.07.02
%
% $Id$


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Process command line arguments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Test number of arguments
if (nargin ~= 2)
    error('Must have two input arguments.');
end

% ---- Check for valid input arguments.
if (size(Fp)~=size(Fc))
    error('Input arguments must be the same size.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert to DP frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Number of detectors.
numberOfDetectors = size(Fp,2);

% ---- Compute rotation needed to reach DP frame.
psi = 1/4*atan(2*(sum(Fp.*Fc,2))./(sum(Fp.*Fp,2)-sum(Fc.*Fc,2)));
psiDP = repmat(psi,[1,numberOfDetectors]);

% ---- Rotate to DP frame.
FpDP = cos(2*psiDP).*Fp + sin(2*psiDP).*Fc;
FcDP = -sin(2*psiDP).*Fp + cos(2*psiDP).*Fc;
Fp = FpDP;
Fc = FcDP;

% ---- Further rotate polarization by pi/4 if |Fp|<|Fc|.
swapIndex = logical(sum(FpDP.^2,2)<sum(FcDP.^2,2));
Fp(swapIndex,:) =  FcDP(swapIndex,:);
Fc(swapIndex,:) = -FpDP(swapIndex,:);
psiDP(swapIndex,:) = psiDP(swapIndex,:) + pi/4;
psi(swapIndex) = psi(swapIndex) + pi/4;

% ---- Done.
return

