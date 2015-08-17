function [xVec, probDens, cumDist] = calcprobdensity(x, xVec)
% 
% CALCPROBDENSITY - Calculate the probability density and cumulative
% distribution function of a sample of random variables.
%
% Usage: [xVec, probDens, cumDist] = calcprobdensity(x, xVec)
%
% x        : Sample
% xVec     : Vector at which the prob. density should be calculated
% probDens : Probability denstity
% cumDist  : Cumulative distribution function
%
% P. Ajith, 24.02.06
%
% $Id: calcprobdensity.m 162 2009-08-14 20:11:29Z ajith $

for i = 1 : length(xVec)
    index = find(x <= xVec(i));
    cumDist(i) = length(index) / length(x);
    
    index = [];
    dx    = (max(xVec) - min(xVec)) / length(xVec);
    index = find(x >= xVec(i) & x < xVec(i) + dx);
    probDens(i) = length(index) / length(x) / dx;    
end