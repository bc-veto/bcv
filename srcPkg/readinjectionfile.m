function [injectionParameters] = readinjectionfile(startTime,blockSize,injectionFileName)
% readinjectionfile - Read simulated signal parameters from a file.
%
%   [injectionParameters] = readinjectionfile(startTime,blockSize,injectionFileName)
%
%  startTime    Scalar.  Start time of the data being analysed.
%  blockSize    Scalar.  Duration (s) of the data being analysed.
%  injectionFileName
%               String.  Name of file specifying simulated signals to be
%               injected.
%
%  injectionParameters
%               Cell array of strings.  Each string contains the parameters
%               for a single GWB or "glitch" simulation.  Only injections
%               in the interval [startTime,startTime+blockSize] are
%               returned.
%
% $Id: readinjectionfile.m 1962 2007-12-18 21:43:18Z jrollins $

%----- Check number of arguments
if (nargin~=3) 
    error('Wrong number of input arguments.')
end

%----- Scan files contents into cell array of strings, one string per line
%      in file.
innerContents = dataread('file', injectionFileName,'%s','delimiter','\n');

% Wrap for compatibility with code
fileContents{1} = innerContents;

%----- Retrieve GPS time of each injection and determine which injections
%      fall in time interval of interest.
nInjection = size(fileContents{1},1);
keepInjection = [];  %-- row numbers of injections to keep
for jInjection = 1:nInjection
    [gps_s rest] = strtok(fileContents{1}{jInjection});
    gps_ns = strtok(rest);
    gps = str2num(gps_s) + 1e-9 * str2num(gps_ns);
    % if ( (gps >= startTime-1) & (gps <= startTime+blockSize+1) )  %-- KLUDGE; actually want all injections overlapping interval.
    if ( (gps >= startTime) & (gps <= startTime+blockSize) )
        keepInjection = [ keepInjection; jInjection];
    end
end

%----- Write desired injections to output cell array.
injectionParameters = {};
for jInjection=1:length(keepInjection)
    injectionParameters{jInjection,1} = fileContents{1}{keepInjection(jInjection)};
end

%----- Done.
return
