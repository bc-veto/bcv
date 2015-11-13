function allyaxis(y1, y2)

% Set all the yaxis ranges on the current
% figure.
%
% usage:  allyaxis(y1, y2)
%
% M Hewitson 17-11-02
%
% $Id: allyaxis.m,v 1.1 2002/11/17 12:48:20 hewitson Exp $
%




c = get(gcf, 'children');

for k=1:length(c)

  t = get(c(k), 'Tag');
  if isempty(t)
    set(c(k), 'YLim', [y1 y2]);
  end
   
end

% 
% END