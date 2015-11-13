function vetoanalysis(frameCache, chanHName, chanXName, frameTypeChanH, ...
    	frameTypeChanX, samplFreqH, samplFreqX, highPassCutoff, TriggerHList, ...
	TriggerXList, couplingModel, transFnXtoH, analysisStartTime, analysisEndTime, ...
	timeShift, outDir, logFid, debugLevel)

% VETOANALYSIS - Perform veto analysis using known instrumental couplings
% on a set of triggers in the GW channel making use of a set of triggers in 
% an instrumental channel X and the transfer function from channel X to the 
% GW channel. 
% 
% Usage: vetoanalysis(TriggerHList, TriggerXList, transFnXtoH, ...
%                     numberOfChannels, timeShift)
%
% TriggerHList - A structure containing the start times, end times, and
%                   central times of the triggers in channel H
% TriggerXList - A structure containing the start times, end times, and
%                   central times of the triggers in channel X
% transFnXtoH     - A structure containing the transfer functions (frequency
%                   and amplitude) from channel X to channel H
% analysisStartTime       - GPS start time lower bound of candidate events
% analysisEndTime         - GPS end time upper bound of candidate events
% timeShift       - The number of seconds by which the coincident triggers will
%                   be time shifted after identification.
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 06-07-09
% P. Ajith <ajith@caltech.edu>, 23-07-09
% 
% $Id: vetoanalysis.m 319 2010-01-05 10:03:13Z isogait $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Find Coincident Triggers Between Channels H And X            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxCoincCentTime = max([length(TriggerHList.centralTime) ...
    length(TriggerXList.centralTime)]);

COINC_TIME_WINDOW = 0.5;  % time-window for identifying coincidences between channels H and X
MAX_LENGTH_DATA_SEG = 64; % maximum lenght of one data segment used for the analysys

% Set the segment length for detecting coincident triggers between channels
% H and X.
segLength = 3600;

% Set the uniqueness argument for detecting coincident triggers between
% channels H and X.
uniqueArgument = 'nonunique';

% Identify coincident triggers between channels H and X using gps central
% times.
fprintf(logFid, ...
	'LOG: Finding coincidences... Num trigs H = %d, Num trigs X = %d\n', ...
	length(TriggerHList.centralTime), length(TriggerXList.centralTime));

[coincTrigH, coincTrigX] = mcoinc(maxCoincCentTime, ...
    TriggerHList.centralTime', TriggerXList.centralTime' + timeShift, ...
    COINC_TIME_WINDOW, segLength, uniqueArgument);

% Get the central times of the coincident triggers in channels H
% and X.
trigHCentTimeVec = TriggerHList.centralTime(coincTrigH)';
trigXCentTimeVec = TriggerXList.centralTime(coincTrigX)';

trigHStartTimeVec = TriggerHList.startTime(coincTrigH)';
trigXStartTimeVec = TriggerXList.startTime(coincTrigX)';

trigHEndTimeVec = TriggerHList.endTime(coincTrigH)';
trigXEndTimeVec = TriggerXList.endTime(coincTrigX)';

trigHCentFreqVec = TriggerHList.centralFrequency(coincTrigH)';
trigXCentFreqVec = TriggerXList.centralFrequency(coincTrigX)';

trigHSignificVec = TriggerHList.triggerSignificance(coincTrigH)';
trigXSignificVec = TriggerXList.triggerSignificance(coincTrigX)';

clear TriggerHList TriggerXList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Read Data Segments Using A Time-Window Around Coincident Triggers    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the sample frequencies for channels H and X are the same. If
% the sample frequencies for channels H and X are not the same, display an
% error message.
if samplFreqH ~= samplFreqX
    error('ERROR: samplFreqH DOES NOT equal samplFreqX');
else
    samplFreq = samplFreqH;
end

% the time-shift applied to the X/Y data. In the case of biliniear coupling
% there are two channels (X and Y); thus a vector of length 2
if strcmp(couplingModel, 'linear')
	timeShiftX  = timeShift;
elseif strcmp(couplingModel, 'bilinear')
	timeShiftX  = [timeShift timeShift];
end

% Check that the number of triggers in channel H that are coincident in
% channel X (coincTrigH) IS the same as the number of triggers in channel X
% that are coincident in channel H (coincTrigX).
if size(coincTrigH) == size(coincTrigX) & length(coincTrigH) > 0
    
    % Read data for timeShift from framecache.
	fprintf(logFid, ...
		'LOG: Performing veto analysis for time shift %d...\n', timeShift);
	fprintf(logFid, 'LOG: Number of coincident triggers: %d...\n', ...
		length(coincTrigH));
    
    % Open a file for writing the cross-correlation statistics for each
    % timeshift veto analysis.
    outFileString = sprintf('/corrstat_timeshift%d_seg%d-%d.dat', ...
        timeShift, analysisStartTime, analysisEndTime);
    outFileName = strcat(outDir, outFileString);
    outFid = fopen(outFileName, 'w+');
    
    for coincIndex = 1 : length(coincTrigH)
        
        % Get the parameters of THIS coincident trigger
        trigHCentTime = trigHCentTimeVec(coincIndex);
        trigXCentTime = trigXCentTimeVec(coincIndex);

        trigHStartTime = trigHStartTimeVec(coincIndex);
        trigXStartTime = trigXStartTimeVec(coincIndex);
        
        trigHEndTime = trigHEndTimeVec(coincIndex);
        trigXEndTime = trigXEndTimeVec(coincIndex);
        
        trigHCentFreq = trigHCentFreqVec(coincIndex);
        trigXCentFreq = trigXCentFreqVec(coincIndex);
        
        trigHSignific = trigHSignificVec(coincIndex);
        trigXSignific = trigXSignificVec(coincIndex);
        
        trigHDuration = trigHEndTime - trigHStartTime;
        trigXDuration = trigXEndTime - trigXStartTime;

	% Find the segment of data used for analysis.
	segStartTime = min([trigHStartTime (trigXStartTime + timeShift)]);
	segEndTime = max([trigHEndTime (trigXEndTime + timeShift)]);
	meanTrigCentTime = (segStartTime + segEndTime) / 2;
	totalDuration = segEndTime - segStartTime;
        
        % Calculate the total number of samples, rounded to the closest
        % integer power of two.
        totalNumberSamples = roundtopowertwo(totalDuration * ...
            samplFreq, 1024);
        
        % Calculate the length of the data segment used for veto analysis.
        segmentDuration = totalNumberSamples / samplFreq;
        
        % Find the start and end times of a small segment of data that we
        % want to analyze.
        segStartTime = meanTrigCentTime - segmentDuration / 2;
        segEndTime = meanTrigCentTime + segmentDuration / 2;
        
        % Round the start and end times of a small segment of data that we
        % want to analyze to an integer gps time.
        % IMPORTANT NOTE: One second of data in the beginning is used to 
        % train the high-pass filter and hence should not be used for the
        % analysis (this explains "floor(segStartTime) - 1" below).
        wreadStartTime = floor(segStartTime) - 1;
        wreadEndTime = ceil(segEndTime) + 1;
        
        if debugLevel >= 2

            % print values of different variables
            printvar('--- Analysis Window ---',[],'analysisStartTime',...
                analysisStartTime, 'analysisEndTime', analysisEndTime, ....
                'wreadStartTime', wreadStartTime, 'wreadEndTime', ....
                wreadEndTime, 'segStartTime', segStartTime, 'segEndTime', ...
                segEndTime, '--- Trigger Parameters ---',[],'trigHStartTime', ...
                trigHStartTime, 'trigXStartTime', trigXStartTime, ...
                'trigHEndTime', trigHEndTime, 'trigXEndTime', trigXEndTime, ...
                'trigHCentTime', trigHCentTime, 'trigXCentTime', trigXCentTime, ...
                'trigHCentFreq', trigHCentFreq, 'trigXCentFreq', trigXCentFreq, ...
                'trigHSignific', trigHSignific, ...
                'trigXSignific', trigXSignific, ....
                '--- Veto Analysis Parameters ---', [], 'meanTrigCentTime', ...
                meanTrigCentTime, 'totalDuration', totalDuration);
        end
        
		% check whether the segments length of the time-seriew data is sensible
		if (wreadEndTime - wreadStartTime) > MAX_LENGTH_DATA_SEG
			fprintf(logFid, ...
				'ERROR: Segment length %f is larger than the allowed max length %f ...\n', ...
					wreadEndTime - wreadStartTime, MAX_LENGTH_DATA_SEG);
		end

        if segStartTime >= analysisStartTime && segEndTime <= analysisEndTime && ...
			(segEndTime-segStartTime) <= MAX_LENGTH_DATA_SEG

            % Read this small segment of data for channel H and channel X.
            [dataH, samplFreqH] = ...
                wreaddata(frameCache, chanHName, frameTypeChanH, ...
                wreadStartTime, wreadEndTime, 0, debugLevel);
            
            [dataX, samplFreqX] = ...
                wreaddata(frameCache, chanXName, frameTypeChanX, ...
                wreadStartTime, wreadEndTime, timeShiftX, debugLevel);
            
            % Check for a read error in the channel H data.
            if ~all(samplFreqH)
                fprintf(logFid, ...
                    'ERROR: Cannot load frame data for channel H...\n');
                fprintf(logFid, ...
                    'ERROR: Channel H - wreadStartTime: %f wreadEndTime %f...\n', ...
                    wreadStartTime, wreadEndTime);
            
            % Check for a read error in the channel X data.
            elseif ~all(samplFreqX)
                fprintf(logFid, ...
                    'ERROR: Cannot load frame data for channel X\n');
                fprintf(logFid, ...
                    'ERROR: Channel X - wreadStartTime: %f wreadEndTime %f...\n', ...
                    wreadStartTime + timeShift, wreadEndTime + timeShift);
            
            else
                
				% Create a time vector for the data.
				clear timeH timeX
                timeH = [wreadStartTime : 1/samplFreqH: wreadEndTime-1/samplFreqH];
				for iX = 1:length(samplFreqX)
                	timeX{iX} = [wreadStartTime-timeShift:1/samplFreqX(iX):wreadEndTime-timeShift-1/samplFreqX(iX)];
				end
                
                % plot time-series data, if necessary
                if debugLevel >= 2
                    
                    close all
                    
                    fh1 = plotdata(timeX{1}-min(timeX{1}), dataX{1}, 10, 'lin', 'lin', ...
                        ['t[sec] since ', num2str(min(timeX{1}))], 'x(t)', ...
                        ['Time-series data:' , strrep(chanXName{1},'_','\_')], ...
                        {'raw data','high-passed'},{'m'});
                        
                    fh2 = plotdata(timeH-min(timeH), dataH, 11, 'lin', 'lin', ...
                        ['t[sec] since ', num2str(min(timeH))], 'h(t)', ...
                        ['Time-series data:' , strrep(chanHName,'_','\_')], ...
                        {'raw data','high-passed'},{'c'});

                    if length(dataX) > 1
                        fh3 = plotdata(timeX{2}-min(timeX{2}), dataX{2}, 12, 'lin', 'lin', ...
                            ['t[sec] since ', num2str(min(timeX{2}))], 'y(t)', ...
                            ['Time-series data:' , strrep(chanXName{2},'_','\_')], ...
                            {'raw data','high-passed'},{'g'});
                    end

                end

                % if the sampling frequency of any of the data streams (read from 
                % frame files) is different from the one specified in the configuration
                % files, resample the data. 
			    if samplFreqH ~= samplFreq 
				    dataH = wresample(dataH, samplFreqH, samplFreq);
                end
				if ~all(samplFreqX == samplFreq)
					dataX = wresample(dataX, samplFreqX, samplFreq);
                end

				% Create a time vector for the data.
				clear timeH timeX
                timeH = [wreadStartTime : 1/samplFreq: wreadEndTime-1/samplFreq];
                timeX = [wreadStartTime-timeShift : 1/samplFreq : wreadEndTime-timeShift-1/samplFreq];
                
                % in the case of bilinear coupling multiply the X and Y channels
                % to form a pseudo channel (which is a bilinear combination of 
                % X and Y)
                if strcmp(couplingModel,'bilinear')

					% read out some parameters describing the Y channel
    				segIdx = find(timeX+timeShift >= segStartTime & timeX+timeShift < segEndTime);

					meanY = mean(dataX{2}(segIdx));
					varY  = var(dataX{2}(segIdx));
					maxY  = max(dataX{2}(segIdx));
					minY  = min(dataX{2}(segIdx));
			
                    dataXtmp = dataX{1}.*dataX{2};
					clear dataX
					dataX{1} = dataXtmp;
					clear dataXtmp
				else 
					meanY = 0;
					varY  = 0;
					maxY  = 0;
					minY  = 0;
                end

                % High-pass the data.
                if highPassCutoff > 0
                	tiling.sampleFrequency = samplFreq;
                	tiling.highPassCutoff = highPassCutoff;
                	[dataH] = whighpass(dataH, tiling);
                	[dataX] = whighpass(dataX, tiling);
				end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % compute the correlation coefficient according to different noise
                % coupling models -- core of the analysis
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                switch couplingModel
                case 'linear'

                    [crossCorrXH, crossCorrMaxXH] = linearcouplingcoeff(chanHName, ...
                        chanXName, dataH, dataX, timeH, timeX, transFnXtoH, segStartTime, ...
                        segEndTime, timeShift, samplFreq, logFid, debugLevel);
                case 'bilinear'

                    [crossCorrXH, crossCorrMaxXH, dataHClean] = bilinearcouplingcoeff(chanHName, ...
                        chanXName, dataH, dataX, timeH, timeX, segStartTime, segEndTime, ...
                        timeShift, samplFreq, logFid, 1, debugLevel);

                otherwise  

                    error('ERROR: Unknown coupling model.');
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Print the cross-correlation statistic to a file
                % with the timeshift information.
                fprintf(outFid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', timeShift, crossCorrXH, ...
					crossCorrMaxXH, trigHCentTime, trigXCentTime, trigHCentFreq, trigXCentFreq, ...
					trigHSignific, trigXSignific, trigHDuration, trigXDuration, meanY, varY, ...
					maxY, minY); 

            end     % END if statement checking whether the data read from framefiles is empty.
        end         % END if statement segStartTime > analysisStartTime ....
    end             % END for loop over coincident triggers 
    
    % Close the file for writing the cross-correlation statistics for each
    % timeshift veto analysis.
    fclose(outFid);
else
    fprintf(logFid, ...
        'WARNING: No coincident triggers found for timeshift = %d...\n', ...
        timeShift);
end
