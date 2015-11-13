function allxaxis(x1, x2)

% Set all the xaxis ranges on the current
% figure.
%
% usage:  allxaxis(x1, x2)
%
% M Hewitson 17-11-02
%
% $Id: allxaxis.m,v 1.1 2002/11/17 12:48:20 hewitson Exp $
%




c = get(gcf, 'children');

for k=1:length(c)

  t = get(c(k), 'Tag');
  if isempty(t)
    set(c(k), 'XLim', [x1 x2]);
  end
   
end

% 
% END