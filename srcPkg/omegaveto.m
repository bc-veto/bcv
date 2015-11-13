function [] =  omegaveto(segmentFile, configurationFile, frameCacheFile, ...
        couplingModel, highPassCutoff, trigSignThreshX, trigSignThreshH, ...
        outDir, logFile, timeShiftMin, timeShiftMax, numTimeShifts, debugLevel)
%
% OMEGAVETO - perform veto analysis using different coupling models between 
% the instrumental channel and the gravitational-wave channel. 
% 
% usage: omegaveto(segmentFile, configurationFile, frameCacheFile, ...
%        couplingModel, highPassCutoff, trigSignThreshX, trigSignThreshH, ...
%        outDir, logFile, timeShiftMin, timeShiftMax, numTimeShifts, debugLevel)
%        
% configurationFile - Path name of channel configuration file
% frameCacheFile 	- Path name of framecache file
% outDir 			- Directory to write results
% logFile			- string specifying the name of the log file 
% debugLevel 		- Verboseness of debug level output
%
% The configuration file is an ASCII text file describing the parameters
% for each channel to be analyzed. The entries should have the following
% syntax:
%
% We will assume that the first channel that is read from the configuration
% file is channel H (the h(t) detector output data vector), followed by a
% X channel and Y channel (the "fast" and "slow" aux. channles). 
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 22-06-09
% P. Ajith <ajith@caltech.edu>, 23-07-09
% Tomoki Isogai <isogait@carleton.edu>

% Verify correct number of input arguments.
if nargin < 13

	fprintf('\nVeto analysis pipeline using instrumental coupling models:\n\n');

	fprintf('usage (binary executable from a unix shell): \n');
	fprintf('\tomegaveto <segmentFile> <configurationFile> <frameCacheFile> \\ \n');
    fprintf('\t\t<couplingModel> <highPassCutoff> <trigSignThreshX> <trigSignThreshH> \\ \n');
    fprintf('\t\t<outDir> <logFile> <timeShiftMin> <timeShiftMax> <numTimeShifts> <debugLevel>\n\n') 
	
	fprintf('usage (from matlab command window): \n');
	fprintf('\tomegaveto(segmentFile, configurationFile, frameCacheFile, ...\n');
	fprintf('\tcouplingModel, highPassCutoff, trigSignThreshX, trigSignThreshH, ...\n');
	fprintf('\toutDir, logFile, timeShiftMin, timeShiftMax, numTimeShifts, debugLevel)\n');
	
	fprintf('\nP. Ajith, Aaron B. Pearlman, Tomoki Isogai 2009-2010\n\n');
      
        segmentFile='/home/bernard.hall/DETCHAR/omega_test/test15HigherThresh/segfiles/H1_1078000000_1078100000_ER5_TEST_1078000000_1078020000_segs.txt';

	configurationFile='/home/bernard.hall/DETCHAR/omega_test/test15HigherThresh/configs/H1_1078000000_1078100000_ER5_TEST_H1_PSL-ISS_PDA_OUT_DQ_1024_4096+H1_SUS-ITMY_L3_OPLEV_PIT_OUT_DQ_1_16.conf';

	frameCacheFile='/home/bernard.hall/DETCHAR/omega_test/test15HigherThresh/cache/H-H1_R_CACHE-1078000000.0-20000.0.lcf';

	couplingModel='bilinear';
	highPassCutoff=32;
	trigSignThreshX=7.0;
	trigSignThreshH=7.0;

	outDir='/home/bernard.hall/DETCHAR/omega_test/test15HigherThresh/results/H1_1078000000_1078100000_ER5_TEST-H1-H1_PSL-ISS_PDA_OUT_DQ_1024_4096+H1_SUS-ITMY_L3_OPLEV_PIT_OUT_DQ_1_16/1078000000_1078020000';

	logFile='log.txt';
	timeShiftMin=-180;
	timeShiftMax=180;
	numTimeShifts=51;
	debugLevel=0;

	%return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Defaults                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default configuration parameters.
defaultTriggerListChH = '/archive/home/apearlma/opt/omega-veto/config/E14/FullE14/H1_LSC-DARM_ERR_64_1024.trg';
defaultTriggerListChX = '/archive/home/apearlma/opt/omega-veto/config/E14/FullE14/H1_LSC-MICH_CTRL_64_1024.trg';
defaultTransferFunctionXtoH = '/archive/home/apearlma/opt/omega-veto/config/E14/FullE14/H1_MICH_TF_DC.mat';
defaultSampleFrequency = 4096;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Create Output Directory                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If outDir not specified, make one based on center time.
if isempty(outDir)
    error('# outDir empty.');
end

% Report status.
fprintf('Creating event directory...\n');
fprintf('OutputDirectory: %s\n', outDir);

% Create spectrogram directory.
unix(['mkdir -p ' outDir]);

% Copy configuration file.
unix(['cp ' configurationFile ' ' outDir '/configuration.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Open An Output Log File For Writing                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logFile = [outDir, '/', logFile];
logFid = fopen(logFile, 'w+');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Process Command Line Arguments                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isstr(trigSignThreshH)
    trigSignThreshH = str2num(trigSignThreshH);
end

if isstr(trigSignThreshX)
    trigSignThreshX = str2num(trigSignThreshX);
end

if isstr(highPassCutoff)
    highPassCutoff = str2num(highPassCutoff);
end

if isstr(timeShiftMin)
    timeShiftMin = str2num(timeShiftMin);
end

if isstr(timeShiftMax)
    timeShiftMax = str2num(timeShiftMax);
end

if isstr(numTimeShifts)
    numTimeShifts = str2num(numTimeShifts);
end

if isstr(debugLevel)
    debugLevel = str2num(debugLevel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Read Segment File                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segList = load(segmentFile);

segStartTimeVec = segList(:,2);
segEndTimeVec = segList(:,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Write Log File Header Information                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Report status.
if debugLevel >= 0
    fprintf(logFid, '## Bilinear-coupling-veto analysis [%d - %d]\n', ...
        min(segStartTimeVec), max(segEndTimeVec));
    fprintf(logFid, '## Created by %s on %s at %s\n', ...
        getenv('USER'), datestr(clock, 29), datestr(clock, 13));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Read Configuration File                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Report status.
if debugLevel >= 0
    fprintf(logFid, 'LOG: Reading configuration file %s...\n', ...
        configurationFile);
end

% Enable or disable warnings.
if debugLevel >= 2
    warning on OMEGAVETO:incompleteConfiguration
else
    warning off OMEGAVETO:incompleteConfiguration
end

% Initialize configuration structure.
configuration = cell(2, 1);
channelNumber = 0;
sectionIndex = [];
sectionName = [];
sectionNumber = 0;
sectionStart = [];

% Open configuration file for reading.
configurationFID = fopen(configurationFile, 'r');

if configurationFID <= 0
    error(sprintf('### Unable to open the configuration file:%s', ...
        configurationFile));
end

% Begin loop over configuration file. Read only the first 3 channels 
% if more is present
while ~feof(configurationFID) && channelNumber <= 3

  % Read one line from configuration file.
  configurationLine = fgetl(configurationFID);

  % Remove any comments.
  commentIndices = findstr(configurationLine, '#');
  if ~isempty(commentIndices),
    configurationLine = configurationLine(1 : (commentIndices(1) - 1));
  end

  % Remove leading and trailing blanks.
  configurationLine = fliplr(deblank(fliplr(deblank(configurationLine))));

  % If empty line, skip to the next line.
  if isempty(configurationLine),
    continue;
  end

  % Check for new section.
  if configurationLine(1) == '[',

    % Locate field separator.
    commaIndex = strfind(configurationLine, ',');

    % If field separator not located, report syntax error.
    if isempty(commaIndex),
      error('Syntax error processing configuration file "%s":\n%s\n', ...
            configurationFile, configurationLine);
    end

    % Select first field separator.
    commaIndex = commaIndex(1);

    % Increment section number.
    sectionNumber = sectionNumber + 1;

    % Extract section index.
    sectionIndex{sectionNumber} = configurationLine(2 : commaIndex - 1);

    % Extract section name.
    sectionName{sectionNumber} = configurationLine((commaIndex + 1) : end - 1);

    % Remove leading and trailing blanks.
    sectionIndex{sectionNumber} = ...
        fliplr(deblank(fliplr(deblank(sectionIndex{sectionNumber}))));
    sectionName{sectionNumber} = ...
        fliplr(deblank(fliplr(deblank(sectionName{sectionNumber}))));

    % Standardize section names.
    sectionIndex{sectionNumber} = strrep(sectionIndex{sectionNumber}, ';', ':');
    sectionName{sectionNumber} = strrep(sectionName{sectionNumber}, ';', ':');

    % sectionIndex{sectionNumber} = strrep(sectionIndex{sectionNumber}, ';', '_');
    % sectionName{sectionNumber} = strrep(sectionName{sectionNumber}, ';', '_');

    % Determine initial visibility.
    switch sectionIndex{sectionNumber},
      case 'Context',
        sectionChecked{sectionNumber} = 'Checked';
        sectionDisplay{sectionNumber} = 'Block';
        
      case 'Gravitational',
        sectionChecked{sectionNumber} = 'Checked';
        sectionDisplay{sectionNumber} = 'Block';
        
      otherwise
        sectionChecked{sectionNumber} = 'Checked';
        sectionDisplay{sectionNumber} = 'Block';
    end

    % Record first channel in section.
    sectionStart(sectionNumber) = channelNumber + 1;

    % Continue to next line.
    continue;

  end

  % Check for beginning of new channel configuration.
  if configurationLine == '{'
      % Increment channel number
      channelNumber = channelNumber + 1;
      
      % Allocate space for the configuration structure to store data for
      % the first two channels read from the configuration file.
      if channelNumber <= 2
          % Initialize configuration parameters
          configuration{channelNumber}.channelName = [];
          configuration{channelNumber}.frameType = [];
          configuration{channelNumber}.sampleFrequency = [];
          configuration{channelNumber}.triggerListChH = [];
          configuration{channelNumber}.triggerListChX = [];
          configuration{channelNumber}.transferFunctionXtoH = [];
      end

      % Continue to next line.
      continue;
  
  end

  % Check for end of existing channel configuration.
  if configurationLine == '}'
      % Validate channel configuration.
      if isempty(configuration{channelNumber}.channelName)
          error('Channel name not specified for channel number %d', ...
              channelNumber);
      end
      
      if isempty(configuration{channelNumber}.frameType)
          error('Frame type not specified for channel number %d', ...
              channelNumber);
      end

      if isempty(configuration{channelNumber}.sampleFrequency)
          warning('OMEGAVETO:incompleteConfiguration', ...
              'Sample frequency not specified for channel number %d', ...
              channelNumber);
          
          configuration{channelNumber}.sampleFrequency = ...
              defaultSampleFrequency;
      end
      
      if channelNumber == 1 && ...
              isempty(configuration{channelNumber}.triggerListChH)
          warning('OMEGAVETO:incompleteConfiguration', ...
              ['Channel H trigger list not specified for channel ' ...
              'number %d'], channelNumber);
          
          configuration{channelNumber}.triggerListChH = ...
              defaultTriggerListChH;
      end
      
      if channelNumber == 2 && ...
              isempty(configuration{channelNumber}.triggerListChX)
          warning('OMEGAVETO:incompleteConfiguration', ...
              ['Channel X trigger list not specified for channel ' ...
              'number %d'], channelNumber);
          
          configuration{channelNumber}.triggerListChH = ...
              defaultTriggerListChX;
      end
      
      if channelNumber == 2 && ...
              isempty(configuration{channelNumber}.transferFunctionXtoH)
          warning('OMEGAVETO:incompleteConfiguration', ...
              ['Channel X transfer function not specified for channel ' ...
              'number %d'], channelNumber);
          
          configuration{channelNumber}.transferFunctionXtoH = ...
              defaultTransferFunctionXtoH;
      end
      
      % Continue to next line.
      continue;
      
  end

  % Locate field separator.
  %name_flag=0;
  
  colonIndex = strfind(configurationLine, ':');  %set index to config colon
  fprintf('configurationLine: %s\n', configurationLine);
  %colonIndex = strfind(configurationLine, '-'); % find first string
  %colonIndex = strfind(configurationLine, '_'); % find indices of this string 

  % If field separator not located, report syntax error.
  if isempty(colonIndex),
    error('Syntax error processing configuration file "%s":\n%s\n', ...
          configurationFile, configurationLine);
  end

  colonIndex_2 = numel(colonIndex); % value of 1 means new convention, a value of 2 means old convention.

  % Parse configuration line.
  colonIndex = colonIndex(1); % use the first instance of the character
  parameterName = configurationLine(1 : colonIndex);
  parameterValue = configurationLine((colonIndex + 1) : end);
  
  fprintf('parameterName (1): %s\n', parameterName);
  fprintf('parameterValue (1): %s\n', parameterValue);  

  parameterName = fliplr(deblank(fliplr(deblank(parameterName))));
  parameterValue = fliplr(deblank(fliplr(deblank(parameterValue))));
  
  channel_flag = strcmp(parameterName,'channelName:');

  if (colonIndex_2==1) && (channel_flag == 1), % check for convention, and whether this is channel info
     
     pv_replace_1 = strtok(parameterValue,'_');
     pv_replace_2 = strcat(pv_replace_1,'_');
     parameterValue = strrep(parameterValue,pv_replace_2,strcat(pv_replace_1,':'));
     % i.e., we are simply transforming any new convention values into old convention to then pass on to the rest of the code. In the future, it will probably be better to make the rest of the code handle the newer convention directly.  But this may work for now.

     %name_flag=1;

     parameterValue=regexprep(parameterValue,'_[0-9]+[0-9]',''); % lose end if new convention
     parameterValue=regexprep(parameterValue,'_[0-9]','')
 
     fprintf('configurationLine (2): %s\n', configurationLine);
  end

  fprintf('parameterName (2): %s\n', parameterName);
  fprintf('parameterValue (2): %s\n', parameterValue);

  %name_flag=0;

  % We assume that the first channel read from the configuration file is
  % channel H and the second channel is channel X. We impose a condition to
  % read only the first two channels in the configuration file.
  % Assign parameters based on name.
  switch parameterName
      case 'channelName:'
          configuration{channelNumber}.channelName = ...
              eval(parameterValue);
          
      case 'frameType:'
          configuration{channelNumber}.frameType = ...
              eval(parameterValue);

      case 'sampleFrequency:'
          configuration{channelNumber}.sampleFrequency = ...
              eval(parameterValue);
          
      case 'triggerListChH:'
          if channelNumber == 1
              configuration{channelNumber}.triggerListChH = ...
                  eval(parameterValue);
          end
          
      case 'triggerListChX:'
          if channelNumber == 2
              configuration{channelNumber}.triggerListChX = ...
                  eval(parameterValue);
          end
          
      case 'transferFunctionXtoH:'
          if channelNumber == 2
              configuration{channelNumber}.transferFunctionXtoH = ...
                  eval(parameterValue);
          end
          
      otherwise
          error(['Unknown configuration parameter ' parameterName]);
  end
  
% End loop over channel configuration file.
end

% Close configuration file.
fclose(configurationFID);

% Number of configured channels.
numberOfChannels = length(configuration);

% Number of sections.
numberOfSections = length(sectionName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Load frame cache file, segment list etc. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load frame file cache.
fprintf(logFid, 'LOG: Reading framecache file %s...\n', frameCacheFile);
frameCache = loadframecache(frameCacheFile);

numTrigsH = 0;
numTrigsX = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Loop Over Different Segments                        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSeg = length(segStartTimeVec);
if nSeg ~= length(segEndTimeVec)
	error('Segment start and end time vectors have diff. length. Check seg. file');
end

for iSeg = 1:nSeg

    analysisStartTime = segStartTimeVec(iSeg);
    analysisEndTime = segEndTimeVec(iSeg);

    fprintf(logFid, '# Processing segment # %d. [%d, %d]...\n', ...
        iSeg, analysisStartTime, analysisEndTime);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 Read Trigger Lists For Channels H and X                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set the channel number index for channels H, X and Y.
    channelHIndex = 1;
    channelXIndex = 2;
    channelYIndex = 3;
    
    % Read the trigger list for channel H.
    triggerListChH = configuration{channelHIndex}.triggerListChH;
    
    % Report status.
    if debugLevel >= 0
        fprintf(logFid, 'LOG: Reading channel H trigger list %s...\n', ...
            triggerListChH);
    end
    
    % Open channel H trigger list file for reading.
    triggerListChHFID = fopen(triggerListChH, 'r');
    
    % Count the number of lines in the channel H trigger list file.
    if triggerListChHFID >= 0 
    
    	trigDataMatrixH = load(triggerListChH);

        [p,q]=size(trigDataMatrixH)
	if p > 100
                %trigDataMatrixH = imresize(trigDataMatrixH, [100 q]);
                warning(sprintf('More than 100 H triggers found.  If this causes excessive resource use, try reducing the segment lengths in the .ini file before running setup.'));
                %[p,q]=size(trigDataMatrixH)
        end

    	gpsTriggerHList.startTime = trigDataMatrixH(:,1);
    	gpsTriggerHList.endTime = trigDataMatrixH(:,2);
    	gpsTriggerHList.centralTime = trigDataMatrixH(:,3);
    	gpsTriggerHList.centralFrequency = trigDataMatrixH(:,4);
        % hacked: using SNR instead of significance
    	gpsTriggerHList.triggerSignificance = sqrt(trigDataMatrixH(:,6) - trigDataMatrixH(:,7));
    
    	clear trigDataMatrixH;
    else
    	error(sprintf('### Unable to open the trigger file for channel H\n(%s)', ...
    		triggerListChH));
    end
    
    % Find gpsStart, gpsEnd, and gpsCentral triggers from channel H that are
    % not between analysisStartTime and analysisEndTime.
    triggerIndexChH = find(gpsTriggerHList.centralTime < analysisStartTime ...
        | gpsTriggerHList.centralTime > analysisEndTime ...
        | gpsTriggerHList.triggerSignificance < trigSignThreshH);
    
    % Remove gpsStart, gpsEnd, and gpsCentral triggers from channel H that are
    % not between analysisStartTime and analysisEndTime. We also removed
    % triggers from channel H that are lower than trigSignThreshH in significance.
    gpsTriggerHList.startTime(triggerIndexChH) = [];
    gpsTriggerHList.endTime(triggerIndexChH) = [];
    gpsTriggerHList.centralTime(triggerIndexChH) = [];
    gpsTriggerHList.centralFrequency(triggerIndexChH) = [];
    gpsTriggerHList.triggerSignificance(triggerIndexChH) = [];
    
    clear triggerIndexChH;
    
    % Read the trigger list for channel X.
    triggerListChX = configuration{channelXIndex}.triggerListChX;
    
    % Report status.
    if debugLevel >= 0
        fprintf(logFid, 'LOG: Reading channel X trigger list %s...\n', ...
            triggerListChX);
    end
    
    % Open channel X trigger list file for reading.
    triggerListChXFID = fopen(triggerListChX, 'r');
    
    % Create a counter for the number of lines in the channel X trigger list
    % file.
    numTriggerXLines = 0;
    
    if triggerListChXFID >= 0
   % issue here--trigDataMatrixX ends up being empty.  Fails on first attempt to access matrix. 

        trigDataMatrixX = load(triggerListChX);
        
        [m,n]=size(trigDataMatrixX)
%        if m > 4000
%		trigDataMatrixX = imresize(trigDataMatrixX, [4000 n]);
%                warning(sprintf('Maximum size exceeded.  Reducing X matrix to 4000 x %s',n));
%                [m,n]=size(trigDataMatrixX)
%        end

        trig_x_flag=0

        if (m>0) && (n>0)
		gpsTriggerXList.startTime = trigDataMatrixX(:,1);
        	gpsTriggerXList.endTime = trigDataMatrixX(:,2);
        	gpsTriggerXList.centralTime = trigDataMatrixX(:,3);
        	gpsTriggerXList.centralFrequency = trigDataMatrixX(:,4);
	        % hacked: using SNR instead of significance
    		gpsTriggerXList.triggerSignificance = sqrt(trigDataMatrixX(:,6) - trigDataMatrixX(:,7));
	        %gpsTriggerXList.triggerSignificance = trigDataMatrixX(:,8);
    
        	clear trigDataMatrixX;
	else
        	% error(sprintf('Matrix is empty!'));
                warning(sprintf('Matrix is empty!...No X triggers...Continuing...'));
		trig_x_flag=1
        end
    else
    	error(sprintf('### Unable to open the trigger file for channel X\n(%s)', ...
                    triggerListChX));
    end
    
    % Find gpsStart, gpsEnd, and gpsCentral triggers from channel X that are
    % not between analysisStartTime and analysisEndTime.
    if trig_x_flag==0
		triggerIndexChX = find(gpsTriggerXList.centralTime < analysisStartTime ...
        	| gpsTriggerXList.centralTime > analysisEndTime ...
        	| gpsTriggerXList.triggerSignificance < trigSignThreshX);
    
    % Remove gpsStart, gpsEnd, and gpsCentral triggers from channel X that are
    % not between analysisStartTime and analysisEndTime. We also removed
    % triggers from channel X that are lower than trigSignThreshX in significance.
    	gpsTriggerXList.startTime(triggerIndexChX) = [];
    	gpsTriggerXList.endTime(triggerIndexChX) = [];
    	gpsTriggerXList.centralTime(triggerIndexChX) = [];
    	gpsTriggerXList.centralFrequency(triggerIndexChX) = [];
    	gpsTriggerXList.triggerSignificance(triggerIndexChX) = [];	

    clear triggerIndexChX;
    
    end

    numTrigsH = numTrigsH + length(gpsTriggerHList.centralTime);
    
    if trig_x_flag==0
   	 numTrigsX = numTrigsX + length(gpsTriggerXList.centralTime);
 
    % read the transfer function in the case of linear coupling
    	if strcmp(couplingModel, 'linear')
    
        	transferFunctionXtoH = configuration{channelXIndex}.transferFunctionXtoH;

        % Read the transfer function for channel X to H.
        	if strcmp (transferFunctionXtoH, 'null') == 0
            		transFnXtoHRead = load(transferFunctionXtoH);
            		transFnXtoH.frequency = transFnXtoHRead.auxff;
            		transFnXtoH.Txh = transFnXtoHRead.trans_func;

        % no transfer function. fill it with ones 
        	else
            		transFnXtoH.frequency = linspace(0.1, 1000, 1001)';
            		transFnXtoH.Txh = ones(size(transFnXtoH.frequency));
        	end
    
    	else
        	transFnXtoH = [];
    	end
    
    % Isolate the channel name for channel H.
    	chanHName = configuration{channelHIndex}.channelName;
    end
    % Isolate the frame type for channel H.
    frameTypeH = configuration{channelHIndex}.frameType;
    
    % Isolate the sample frequency for channel H.
    samplFreqH = configuration{channelHIndex}.sampleFrequency;
    
    if trig_x_flag==0
    % Isolate the parameters for channel X.
    	chanXName{1}  = configuration{channelXIndex}.channelName;
    	frameTypeX{1} = configuration{channelXIndex}.frameType;
    	samplFreqX(1) = configuration{channelXIndex}.sampleFrequency;
    end
    % Isolate the parameters for channel Y.
    if strcmp(couplingModel, 'bilinear')
    	chanXName{2} = configuration{channelYIndex}.channelName;
    	frameTypeX{2} = configuration{channelYIndex}.frameType;
    	samplFreqX(2) = configuration{channelYIndex}.sampleFrequency;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Perform Timeshift And Veto Analysis                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % construct a time-shift vector, make sure that the zero lag is inclued 
    % perform the veto analysis for each time shift (incl. zer lag)

    if trig_x_flag==0
    	timeShiftVec = linspace(timeShiftMin, timeShiftMax, numTimeShifts);
    	timeShiftVec = unique([0 round(timeShiftVec)]);

    	for iTimeShift = 1:length(timeShiftVec)

        	timeShift = timeShiftVec(iTimeShift);
                try
    			vetoanalysis(frameCache, chanHName, chanXName, frameTypeH, ...
            			frameTypeX, samplFreqH, samplFreqX, highPassCutoff, ...
    				gpsTriggerHList, gpsTriggerXList, couplingModel, transFnXtoH, ...
    				analysisStartTime, analysisEndTime, timeShift, outDir, logFid, ...
    			debugLevel);
		catch
			warning(sprintf('Time Shift %s failed!',timeShift));	
		end
    	end

    % clear some memory
    	clear gpsTriggerHList gpsTriggerXList;
    else
        clear gpsTriggerHList;
    end
end %% END loop over segments %%

    %end

if trig_x_flag==1
        numTrigsX=-1;
        %clear gpsTriggerHList;
end

% Open text summary file.
textSummaryFID = fopen([outDir '/summary.txt'], 'w');

% print a summary of the analysis
fprintf(textSummaryFID, '# Summary file of the veto analysis\n');
fprintf(textSummaryFID, 'createdBy : %s\n', getenv('USER'));
fprintf(textSummaryFID, 'createdOn : %s\n', datestr(clock));
fprintf(textSummaryFID, 'analysisStartTime UTC : %s\n', gps2utc(analysisStartTime));
fprintf(textSummaryFID, 'analysisEndTime UTC : %s\n', gps2utc(analysisEndTime));
fprintf(textSummaryFID, 'configurationFile : %s\n', configurationFile);
fprintf(textSummaryFID, 'frameCacheFile : %s\n', frameCacheFile);
fprintf(textSummaryFID, 'couplingModel : %s\n', couplingModel);
fprintf(textSummaryFID, 'highPassCutoff : %3.2f\n', highPassCutoff);
fprintf(textSummaryFID, 'outDir : %s\n', outDir);
fprintf(textSummaryFID, 'logFile : %s\n', logFile);
fprintf(textSummaryFID, 'debugLevel : %d\n', debugLevel);
fprintf(textSummaryFID, 'numTrigsH: %d\n', numTrigsH); 
fprintf(textSummaryFID, 'numTrigsX: %d\n', numTrigsX);

% Close text summary file.
fclose(textSummaryFID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Exit                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Report completion.
if debugLevel >= 0
    fprintf(logFid, 'Finished on %s at %s\n', datestr(clock, 29), ...
        datestr(clock, 13));
end

% Close all figures.
close all;

% Close the output log file.
fclose(logFid);

% Return to calling function.
return;
