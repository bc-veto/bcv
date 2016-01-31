function yaxis(y1,y2)


% Set the Y axis range of the current figure
% 
% 
% usage:  yaxis(x1,x2)
% 
% 
% M Hewitson  09-10-01
% 
% $Id: yaxis.m,v 1.1 2001/11/21 08:56:44 martin Exp $
% 
% 

set(gca, 'YLim', [y1 y2]);

% END
