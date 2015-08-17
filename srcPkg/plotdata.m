function [fh] = plotdata(xData, yData, figNum, xScale, yScale, ...
                xLabel, yLabel, titlStr, legVec, colVec, subPlot)
% 
% PLOTDATA is a general-purpose plotting function which will make 2D plots 
% from vectors or cell arrays. 
% 
% plotdata(xData, yData, figNum, xScale, yScale, xLabel, yLabel, ...
%    titlStr, legVec, colVec)
% 
% P. Ajith, 22 Sep 2009
% 
% $Id:$




fh = figure(figNum);

if nargin > 10
    subplot(subPlot)
end
hold on

% X and Y are cell arrays
if iscell(xData) & iscell(yData)
    if length(xData) == length(yData) & length(xData) <= length(colVec)
        for i=1:length(xData)
            plot(xData{i}, yData{i}, colVec{i})
        end
    else
        error('### unequal length for the cell arrays xData, yData or colVec');
    end

% X is a normal vector and Y a cell array
elseif ~iscell(xData) & iscell(yData)
    if length(yData) <= length(colVec)
        for i=1:length(yData)
            plot(xData, yData{i}, colVec{i})
        end
    else
        error('### colVec does not have enough elements');
    end

% X and Y are normal vectors
elseif ~iscell(xData) & ~iscell(yData)
    plot(xData, yData, colVec{1})

% mixed stuff
else
    error('### either xData and yData should be cell arrays, or both should be vectors');
end

% grids, labels etc
grid on
xscale(xScale)
yscale(yScale)
xlabel(xLabel)
ylabel(yLabel)
title(titlStr)
legend(legVec)

