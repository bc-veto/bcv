function [] = printvar(varargin)
% 
% PRINTVAR - print the values of variables on screen. 
% 
% usage: printvar('variableName1', value, 'variableName2', value ...)
% 
% P. Ajith. 29 Sep 2009
% 
% $Id:$


if nargin >= 2
    for i=1:2:nargin
        if ischar(varargin{i})
            if isnumeric(varargin{i+1})
                fprintf('\t %s:\t%s\n', varargin{i}, num2str(varargin{i+1}));
            else
                error('Not a number');
            end
        else
            error('Not a string');
        end
    end
end
