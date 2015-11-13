function wwriteevents(events, fileNames, fieldNames, formatType)
% WWRITEEVENTS Write significant tiles properties to a trigger files
%
% WWRITEEVENTS writes the specified properties of significant tiles to the
% specified trigger files.  Subsequent calls to WWRITEEVENTS can be used to
% append new tile properties to an existing files.  WWRITEEVENTS can write
% results in both TXT and XML formats.
%
% usage: wwriteevents(events, fileNames, fieldNames, formatType);
%
%   events              cell array of significant tile properties
%   fileNames           cell array of file names to append events to
%   fieldNames          cell array of strings of field names to include
%   formatType          format type ('txt' or 'xml') of trigger file
%
% WWRITEEVENTS expects a cell array of event structures with one cell per
% channel as returned by wthreshold or wselect.  The events structure must
% contain the following fields, which describe the properties of statistically
% significant tiles.  Additional optional fields such as amplitude, phase, or
% coherent transform properties may also be specified.
%
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   duration             duration of tile [seconds]
%   bandwidth            bandwidth of tile [Hz]
%   normalizedEnergy     normalized energy of tile []
%
% The field names argument is a cell array of property field names to write
% to the file.  A field name that does not exist in the trigger structure
% will be recorded as zeros in the resulting file.
%
% See also WSEARCH, WEVENT, WTRANSFORM, WTHREHSOLD, and WSELECT.

% Leo C. Stein <lstein@ligo.caltech.edu>
% Shourov K. Chatterji <shourov@ligo.caltech.edu>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 4, nargin));

% apply default arguments
if (nargin < 3) || isempty(fieldNames),
  fieldNames = {'time', 'frequency', 'duration', 'bandwidth', ...
                'normalizedEnergy'};
end
if (nargin < 4) || isempty(formatType),
  formatType = 'unspecified';
end

% force cell arrays
events = wmat2cell(events);
fileNames = wmat2cell(fileNames);

% force one dimensional cell arrays
events = events(:);
fileNames = fileNames(:);

% determine number of channels
numberOfChannels = length(events);

% determine number of field names
numberOfFields = length(fieldNames);

% force lower case format type
formatType = lower(formatType);

% default format type
if strcmp(formatType, 'unspecified'),
  if strcmp(fileNames{1}(end - 3 : end), '.xml')
    formatType = 'xml';
  else
    formatType = 'txt';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate number of file names
if length(fileNames) ~= numberOfChannels,
  error('number of fileNames is inconsistent with number of channels');
end

% validate events structures
for channelNumber = 1 : numberOfChannels,
  if ~isfield(events{channelNumber},'id') || ...
          ~strcmp(events{channelNumber}.id, ...
                  'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

% validate format type
if ~ismember(formatType, {'txt', 'xml'}),
  error(['unknown format type "' formatType '"']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             file format details                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% string for invalid fields
invalidFieldString = '0';

% define field separator character
switch formatType,
  case 'txt',
    fieldDelimiter = ' ';
  case 'xml',
    fieldDelimiter = ',';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                 open file                                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % extract file name
  fileName = fileNames{channelNumber};

  % check if the file already exists
  writeHeaderFlag = ~exist(fileName, 'file');

  % open file for appending
  fileDescriptor = fopen(fileName, 'a');

  % check for error opening file
  if fileDescriptor == -1,
    error(['cannot open file "' fileName '" for appending']);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          write header information                          %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % if file is newly created
  if writeHeaderFlag,

    % switch on format type
    switch formatType

      % handle txt format
      case 'txt',

        % write field names
        for fieldNumber = 1 : numberOfFields,
          fprintf(fileDescriptor,'%% %s\n', fieldNames{fieldNumber});
        end
     
      % handle xml format
      case 'xml',
      
        % write beginning of header
        fprintf(fileDescriptor, ...
                '<?xml version=''1.0'' encoding=''utf-8'' ?>\n');
        fprintf(fileDescriptor, ...
                '<LIGO_LW>\n');
        fprintf(fileDescriptor, ...
                '<Table Name="wsearch:events:table">\n');

        % write field names
        for fieldNumber = 1 : numberOfFields,
          fprintf(fileDescriptor, ...
                  '<Column Name="wsearch:%s" Type="real_8"/>\n', ...
                  fieldNames{fieldNumber});
        end

        % write end of header
        fprintf(fileDescriptor, ...
                '<Stream Name="wsearch:table" Type="Local" Delimiter=",">\n');

    % end switch of format type
    end
    
  % end test for newly created file
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         construct format string                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % initialize format string
  formatString = '';

  % begin loop over field names
  for fieldNumber = 1 : numberOfFields,

    % extract field name
    fieldName = fieldNames{fieldNumber};

    % test for availability of field
    if ~isfield(events{channelNumber}, fieldName),
      formatString = [formatString invalidFieldString fieldDelimiter];
      continue;
    end
    
    % set field format based on field name
    switch fieldName,
      case {'time', 'blockStartTime', 'blockStopTime'},
        formatString = [formatString '%#020.9f' fieldDelimiter];
      case {'livetime'},
        formatString = [formatString '%#.9f' fieldDelimiter];
      case {'probSignal', 'probGlitch'},
        formatString = [formatString '%.16g' fieldDelimiter];  
      otherwise,
        if ischar(events{channelNumber}.(fieldName)),
          formatString = [formatString '%s' fieldDelimiter];
        else
          formatString = [formatString '%#09.3e' fieldDelimiter];
        end
    end

  % end loop over field names
  end

  % remove trailing delimiter and end newline
  formatString = [formatString(1 : end - 1) '\n'];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                            write events to file                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % determine number of triggers
  if ischar(events{channelNumber}.(fieldNames{1})),
    numberOfEvents = 1;
  else
    numberOfEvents = length(events{channelNumber}.(fieldNames{1}));
  end

  % begin loop over events
  for eventNumber = 1 : numberOfEvents,

    % begin loop over field names
    for fieldNumber = 1 : numberOfFields,

      % extract field name
      fieldName = fieldNames{fieldNumber};

      % test for field availability
      if ~isfield(events{channelNumber}, fieldName),
          continue;
      end

      % make the data array
      if ischar(events{channelNumber}.(fieldName)),
        data{fieldNumber} = events{channelNumber}.(fieldName);
      elseif iscell(events{channelNumber}.(fieldName)),
        data{fieldNumber} = events{channelNumber}.(fieldName){eventNumber};
      else
        data{fieldNumber} = events{channelNumber}.(fieldName)(eventNumber);
      end

    % end loop over field names
    end

    % write data to file
    fprintf(fileDescriptor, formatString, data{:});

  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                close file                                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % close file
  fclose(fileDescriptor);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end
