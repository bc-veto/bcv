function xscale(scale)


% Set the X scale of the current axis
% 
% 
% usage:  xscale('scale')
% 
% scale = 'lin' or 'log';
% 
% M Hewitson  02-05-04
% 
% $Id: xaxis.m,v 1.1 2001/11/21 08:56:44 martin Exp $
% 
% 

set(gca, 'XScale', scale);

% END