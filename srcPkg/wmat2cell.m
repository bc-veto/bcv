function m = wmat2cell(m, extratest)
% WMAT2CELL Wrapper for commonly encountered form of mat2cell
%
% The following two calls are completely equivalent:
%
%   if extratest && ~iscell(m),
%     m =  mat2cell(m, size(m, 1), size(m, 2));
%   end
%
%   m = wmat2cell(m, extratest);
%
% If extratest is not specified, it is taken to be true.
%
% See also MAT2CELL.

% Leo C. Stein
% lstein@ligo.caltech.edu

% $Id: wmat2cell.m 2263 2009-08-24 19:54:40Z jrollins $

error(nargchk(1,2,nargin));

if nargin < 2,
  extratest = true;
end

if extratest && ~iscell(m),
  m = mat2cell(m, size(m,1), size(m,2));
end
