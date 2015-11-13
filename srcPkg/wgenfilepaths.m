function outputFiles = wgenfilepaths(channelNames, timeString, ...
                                     outputDirectory, format, outputTypes)
% WGENFILEPATHS create trigger file paths to conform to the LIGO lightweight
% naming convention:
%
%   https://www.lsc-group.phys.uwm.edu/daswg/docs/technical/T050017.txt

% verify correct number of input arguments
error(nargchk(3, 5, nargin));

if (nargin < 4) || isempty(format),
  format = 'txt';
end
if (nargin < 5) || isempty(outputTypes),
  outputTypes = {'DOWNSELECT', 'CLUSTER', 'COINCIDE'};
end

% sort ifo names alphabetically
ifos = cell(1,length(channelNames));
for channelNumber = 1 : length(channelNames),
  ifos{channelNumber} = channelNames{channelNumber}(1:2);
end
% determine network string
ifos = sort(ifos);
networkString = [];
for channelNumber = 1 : length(channelNames),
  networkString = sprintf('%s%s', networkString, ifos{channelNumber});
end

% generate trigger file names
for channelNumber = 1 : length(channelNames),
  for typeNumber = 1 : length(outputTypes),
    outputFiles.(outputTypes{typeNumber}){channelNumber} = ...
        sprintf(...
            '%s/%s-OMEGA_TRIGGERS_%s-%s.%s', ...
            outputDirectory, ...
            channelNames{channelNumber}(1:2), ...
            outputTypes{typeNumber}, ...
            timeString, format);
  end
end

% events file
outputFiles.EVENTS = sprintf(...
    '%s/%s-OMEGA_EVENTS-%s.%s', ...
    outputDirectory, networkString, timeString, format);

% skymap file base name ('@TIME@' to be swapped for event time)
outputFiles.SKYMAP_BASE = sprintf(...
    '%s/%s-OMEGA_SKYMAP-@TIME@.%s', ...
    outputDirectory, networkString, 'txt');
