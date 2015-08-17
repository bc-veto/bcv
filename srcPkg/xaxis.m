function xaxis(x1,x2)


% Set the X axis range of the current figure
% 
% 
% usage:  xaxis(x1,x2)
% 
% 
% M Hewitson  09-10-01
% 
% $Id: xaxis.m,v 1.1 2001/11/21 08:56:44 martin Exp $
% 
% 

set(gca, 'XLim', [x1 x2]);

% END