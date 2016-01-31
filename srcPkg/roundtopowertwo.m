% ROUNDTOPOWERTWO rounds the number to the integer power of 2 closest to
% and greater than N.
%
% Usage: [numRound] = roundtopowertwo(N, numRoundMin)
%
% N - Input number
% numRound - Integer power of 2 closest to and greater than N
% numRoundMin - Minimum desired value of numRound
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 02-07-09
% P. Ajith, 13-07-06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Vetoes for Binary-Coalescence Triggers Using Known Instrumental     %
%                                Couplings                                %
%                                                                         %
%                      Modified By: Aaron B. Pearlman                     %
%                        Mentor: Ajith Parameswaran                       %
%                        Creation Date: 22-06-2009                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [numRound] = roundtopowertwo(N, numRoundMin)

% Calculate an integer power of 2 closest to N.
numRound = 2^ceil(log2(N));

% Calculate an integer power of 2 closest to and greater than N.
if nargin >= 2
    numRoundMin = 2^ceil(log2(numRoundMin));
    
    if numRound < numRoundMin
        numRound = numRoundMin;
    end
end