function x = sosfiltfilt(sos, x)
% SOSFILTFILT Zero-phase filtering with second-order sections
%
% SOSFILTFILT implements zero-phase forward and reverse digital filtering
% using the second-order sections filter implementation.  The resulting
% filter has zero phase distortion, but the squared magnitude of the original
% filter.
%
% usage: y = sosfiltfilt(sos, x);
%
%   sos   second-order sections model for filter
%   x     signal to be filter
%
%   y     resulting filtered signal
%
% If x is a matrix, SOSFILTFILT will filter along the columns of x.
%
% See also SOSFILT and FILTFILT.

% Shourov K. Chatterji
% shourov@ligo.mit.edu

% $Id: sosfiltfilt.m 600 2008-04-04 22:18:21Z jrollins $

% verify correct number of input arguments
error(nargchk(2, 2, nargin));

% determine number of rows
rows = size(x, 1);

% force row vector to column vector
if rows == 1,
  x = x(:);
end

% forward filter data
x = sosfilt(sos, x);

% reverse result
x = flipud(x);

% reverse filter data
x = sosfilt(sos, x);

% reverse result
x = flipud(x);

% return same vector orientation
if rows == 1,
  x = x(:).';
end
