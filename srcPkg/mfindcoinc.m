% MFINDCOINC is a wrapper function for MFINDCOINC.MEX64. MFINDCOINC.MEX64
% is a mex file that looks for coincident events using N parameters.
%
% Usage: function [coincA, coincB] = mfindcoinc(maxcoinc, A, B, win);
%
% Inputs:
%
% maxcoinc - Maximum number of coincidences to return
% A - Matrix of events [m x na]
% B - Matrix of events [m x nb]
% P - Vector of time-windows [m x 1]
%
% The input matrices should have rows of parameters:
%
% e.g.:
%
% >>  A = [0   1.1  4.3   21   43
%          3   5    7.4   8.1  2
%          1   2    3     4    5]
% 
% >>  B = [0.1  2.1 4.2 43
%          3.1  4   2   2
%          1    3   2   5 ]
%
% >> [c1, c2] = mfindcoinc(10, A, B, [0.2; 0.2; 0.1])
%
% Outputs:
%
% c1 - Indices in A that are coincident with an event in B
% c2 - Indices in B that are coincident with an event in A
%
% It is recommended to run 'unique' on the outputs:
%
% >> c1 = unique(c1);
%
% M. Hewitson, 21-06-06
% Aaron B. Pearlman <aaronp1@umbc.edu>, 07-07-09
% P. Ajith, 30-06-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Vetoes for Binary-Coalescence Triggers Using Known Instrumental     %
%                                Couplings                                %
%                                                                         %
%                      Modified By: Aaron B. Pearlman                     %
%                        Mentor: Ajith Parameswaran                       %
%                        Creation Date: 22-06-2009                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%