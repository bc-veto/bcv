% MCOINC is a wrapper function for MFINDCOINC.MEX file. This function assumes
% that the first row of the matrices A and B is a time vector. The matrices
% are split into submatrices of length seglen(secs) for processing. This
% should result in a large speed increase for coincident analysis.
%
% Usage: function [aidx, bidx] = mcoinc(maxcoinc, A, B, P, seglen)
%        function [aidx, bidx] = mcoinc(maxcoinc, A, B, P, seglen, ...
%                                       'nonunique')
%
% This function returns a unique list of indices in A and B. If the
% 'nonunique' flag is given, the indices may not be unique.
%
% A - Parameter matrix of A triggers [m x n]
% B - Parameter matrix of B triggers [m x p]
% P - Coincidence time-window vector [m x 1]
%
% An example usage of this function would be:
%
% [xldx, hldx] = mcoinc(maxCoinc, xTimeVec, hTimeVec, timeWind, 3600, ...
%                       'nonunique')
% 
% In this example, xTimeVec(xldx) will give you the triggers in X that are
% coincident with the triggers in H. hTimeVec(hldx) will give you the
% triggers in H that are coincident with the triggers in X.
%
% timeWind is the time-window (say 0.5 seconds).
%
% maxCoinc is the maximum number of coincidences that you allow (just for
% logistical reasons, like the memory of your computer). If you don't have
% memory limitations, set it to:
%
% maxCoinc = max([length(xTimeVec) length(hTimeVec)]);
%
% M. Hewitson, 27-06-06
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

function [c1, c2] = mcoinc(maxCoinc, A, B, P, seglen, varargin)

% Find the minimum time contained in the first row of matrices A and B.
t0 = floor(min([A(1, :) B(1, :)]));

% Find the maximum time contained in the first row of matrices A and B.
tmax = round(max([A(1, :) B(1, :)]));

% Initialize the output row vectors that will indicate coincident trigger
% columns between matricies A and B.
c1 = [];
c2 = [];

% Set the nonunique flag to false. This will default to the 'unique'
% argument for the varargin input parameter.
nu = false;

% Check the varargin input parameter for the 'nonunique', 'unique', or an
% invalid argument. If varargin is invalid, display an error message.
if nargin > 5
    if strcmp(char(varargin{1}), 'nonunique')
        nu = true;
    elseif strcmp(char(varargin{1}), 'unique')
        nu = false;
    else
        error('Invalid option for varargin in [c1, c2] = mcoinc(maxCoinc, A, B, P, seglen, varargin)');
    end
end

% Split the matrices into submatrices of length seglen(secs) for
% processing. This should result in a large speed increase for coincident
% analysis.
for ts = t0 : seglen : tmax
    
    % Get the index for this segment of the A matrix.
    idxa = find(A(1, :) >= ts & A(1, :) < ts + seglen);
    ATemp = A(:, idxa);
    
    % Get the index for this segment of the B matrix.
    idxb = find(B(1, :) >= ts & B(1, :) < ts + seglen);
    BTemp = B(:, idxb);
    
    % Use the mfindcoinc function to identify coincident trigger columns
    % between matrices ATemp and BTemp.
    [C1, C2] = mfindcoinc(maxCoinc, ATemp, BTemp, P);
    
    % If C1 is not empty, adjust the triggers in A that are coincident with
    % the triggers in B based upon the input parameter varargin,
    % 'nonunique' or 'unique'.
    if ~isempty(C1)
        if nu
            c1  = [c1 C1 + idxa(1) - 1];
        else
            c1  = [c1 unique(C1) + idxa(1) - 1];
        end
    end
    
    % If C2 is not empty, adjust the triggers in B that are coincident with
    % the triggers in A based upon the input parameter varargin,
    % 'nonunique' or 'unique'.
    if ~isempty(C2)
        if nu
            c2 =  [c2 C2 + idxb(1) - 1];
        else
            c2 =  [c2 unique(C2) + idxb(1) - 1];
        end
    end
end