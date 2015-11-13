% INTERPOLATETRANSFN interpolates the supplied transfer function to the
% required frequency resolution.
%
% Usage: [tfFreqInt, tfMagIntp, tfPhaseIntp] = ...
%            interpolatetransfn(tfFreq, tfMag, tfPhase, reqFreqRes)
%
% tfFreq - Frequency vector of the transfer function
% tfMag - Magnitude vector of the transfer function
% tfPhase - Phase vector of the transfer function
% reqFreqRes - Required frequency resolution
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

function [tfFreqIntp, tfMagIntp, tfPhaseIntp] = interpolatetransfn(...
    tfFreq, tfMag, tfPhase, reqFreqRes)

% Compute the required frequency resolution.
tfFRes = tfFreq(2) - tfFreq(1);

% Create a new frequency vector for interpolating the transfer function.
newFreqVec = [min(tfFreq) : reqFreqRes : max(tfFreq)]';

% Use the resample command for interpolating the transfer function.
% tfFreqIntp  = newFreqVec;
% tfMagIntp   = resample(tfMag,tfFRes,reqFreqRes);
% tfPhaseIntp = resample(tfPhase,tfFRes,reqFreqRes);

% Use spline interpolation to interpolate the frequency, magnitude, and
% phase of the transfer function.
tfFreqIntp = interp1(tfFreq, tfFreq, newFreqVec, 'spline');
tfMagIntp = interp1(tfFreq, tfMag, newFreqVec, 'spline');
tfPhaseIntp = interp1(tfFreq, tfPhase, newFreqVec, 'spline');