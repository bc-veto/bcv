function subtractnoise(dbFile, frameCacheFile, chanHName, chanXName, frameTypeChanH, ...
    	frameTypeChanX, samplFreqH, samplFreqX, highPassCutoff,  ...
		couplingModel, transFnXtoH, analysisStartTime, analysisEndTime, ...
		corrThreshMin, corrThreshMax, outDir, logFile, outFile, debugLevel)
%
% SUBTRACTNOISE - Subtract noise transients in the GW channels by projecting
% out certain bilinear combinations of instrumental channels. 
% 
% subtractnoise(dbFile, frameCacheFile, chanHName, chanXName, frameTypeChanH, ...
%	frameTypeChanX, samplFreqH, samplFreqX, highPassCutoff,  ...
%	couplingModel, transFnXtoH, analysisStartTime, analysisEndTime, ...
%	corrThreshMin, corrThreshMax, outDir, logFile, outFile, debugLevel)
%
% P. Ajith <ajith@caltech.edu>, 23-07-09
% 
% $Id: subtractnoise.m 319 2010-01-05 10:03:13Z isogait $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Find Coincident Triggers Between Channels H And X            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% required time shift 
timeShift = 0;

% make output directiries 
system(['mkdir -p ' outDir]);
system(['mkdir -p ' outDir '/plots']);
logFid = fopen([outDir '/' logFile], 'w+');
outFile = [outDir '/' outFile];

params = 'tau, r, rMax, trigHCentTime, trigXCentTime, trigHDuration, trigXDuration';

% get data from the sqlite database
mksqlite('open', dbFile);
Data = mksqlite(['select ' params ' from data where tau == 0 and abs(r) >= ' num2str(corrThreshMin) ' and abs(r) <= ' num2str(corrThreshMax)]);
mksqlite('close');

timeShiftVec 	  = [Data.tau];
crossCorrXHVec	  = [Data.r];
crossCorrMaxXHVec = [Data.rMax];
trigHCentTimeVec  = [Data.trigHCentTime];
trigXCentTimeVec  = [Data.trigXCentTime];
trigHDurationVec  = [Data.trigHDuration];
trigXDurationVec  = [Data.trigXDuration];

% Check that the sample frequencies for channels H and X are the same. If
% the sample frequencies for channels H and X are not the same, display an
% error message.
if samplFreqH ~= samplFreqX
    error('ERROR: samplFreqH DOES NOT equal samplFreqX');
else
    samplFreq = samplFreqH;
end

% open a mat file to print results 
[stat, host] = system('hostname');
[stat, user] = system('whoami');
Info.date = datestr(now);
Info.machine = cellstr(host);
Info.createdBy = cellstr(user);
Info.codeVersion = '$Id: $';

% read frame cache file
% Load frame file cache.
fprintf(logFid, 'LOG: Reading framecache file %s...\n', frameCacheFile);
frameCache = loadframecache(frameCacheFile);

% find triggers in the interval [analysisStartTime, analysisEndTime]
intervalIdx = find(trigHCentTimeVec >= analysisStartTime & trigHCentTimeVec ...
	<= analysisEndTime & timeShiftVec == timeShift); 

if length(intervalIdx) >= 1

	crossCorrXHVec	  = crossCorrXHVec(intervalIdx);
	crossCorrMaxXHVec = crossCorrMaxXHVec(intervalIdx);
	trigHCentTimeVec  = trigHCentTimeVec(intervalIdx);
	trigXCentTimeVec  = trigXCentTimeVec(intervalIdx);
	trigHDurationVec  = trigHDurationVec(intervalIdx);
	trigXDurationVec  = trigXDurationVec(intervalIdx);

    fprintf(logFid, ...
        ['LOG: Total num. triggers ' num2str(length(crossCorrXHVec)) ' ... \n']); 

	% Check that the sample frequencies for channels H and X are the same. If
    for coincIndex = 1 : length(crossCorrXHVec);
        
    	fprintf(logFid, ...
        	['LOG: Processing trigger ' num2str(coincIndex) ' ... \n']); 

        % Get the parameters of THIS coincident trigger
        trigHCentTime = trigHCentTimeVec(coincIndex);
        trigXCentTime = trigXCentTimeVec(coincIndex);

        trigHDuration = trigHDurationVec(coincIndex);
        trigXDuration = trigXDurationVec(coincIndex);

        trigHStartTime = trigHCentTime-trigHDuration/2;
        trigXStartTime = trigXCentTime-trigXDuration/2;
        
        trigHEndTime = trigHCentTime+trigHDuration/2;
        trigXEndTime = trigXCentTime+trigXDuration/2;
        
        % Find the mean of the central times.
        meanTrigCentTime = (trigHCentTime + trigXCentTime) / 2;
        
        % Calculate the total duration of the two triggers.
        totalDuration = max([trigHEndTime (trigXEndTime)]) ...
            - min([trigHStartTime (trigXStartTime)]);
        
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
        wreadStartTime = floor(segStartTime) - 2;
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
                '--- Veto Analysis Parameters ---', [], 'meanTrigCentTime', ...
                meanTrigCentTime, 'totalDuration', totalDuration, 'crossCorrXH', ...
				crossCorrXHVec(coincIndex));
        end
        
        if segStartTime >= analysisStartTime && ...
                segEndTime <= analysisEndTime
            
            % Read this small segment of data for channel H and channel X.
            [dataH, samplFreqH] = ...
                wreaddata(frameCache, chanHName, frameTypeChanH, ...
                wreadStartTime, wreadEndTime, 0, debugLevel);
            
            [dataX, samplFreqX] = ...
                wreaddata(frameCache, chanXName, frameTypeChanX, ...
                wreadStartTime, wreadEndTime, [0, 0], debugLevel);

			Data(coincIndex).dataH = dataH;
			Data(coincIndex).dataX = dataX;
            
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
                    wreadStartTime, wreadEndTime);
            
            else
                
				% Create a time vector for the data.
				clear timeH timeX
                timeH = [wreadStartTime : 1/samplFreqH: wreadEndTime-1/samplFreqH];
				for iX = 1:length(samplFreqX)
                	timeX{iX} = [wreadStartTime:1/samplFreqX(iX):wreadEndTime-1/samplFreqX(iX)];
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

				Data(coincIndex).dataH_resampl = dataH;
				Data(coincIndex).dataX_resampl = dataX;

                % in the case of bilinear coupling multiply the X and Y channels
                % to form a pseudo channel (which is a bilinear combination of 
                % X and Y)
                if strcmp(couplingModel,'bilinear')
                    dataXtmp = dataX{1}.*dataX{2};
					clear dataX
					dataX{1} = dataXtmp;
					clear dataXtmp
                end
                
                % High-pass the data.
                if highPassCutoff > 0
                	tiling.sampleFrequency = samplFreq;
                	tiling.highPassCutoff = highPassCutoff;
                	[dataH] = whighpass(dataH, tiling);
                	[dataX] = whighpass(dataX, tiling);
				end

				Data(coincIndex).dataH_resampl_hp = dataH;
				Data(coincIndex).dataX_resampl_hp = dataX;

				% Create a time vector for the data.
				clear timeH timeX
                timeH = [wreadStartTime: 1/samplFreq: wreadEndTime-1/samplFreq];
                timeX = [wreadStartTime: 1/samplFreq: wreadEndTime-1/samplFreq];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				% project out the noise from the pseudo instrumental channels 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                switch couplingModel
                case 'linear'

                    [crossCorrXHNew, crossCorrMaxXHNew] = linearcouplingcoeff(chanHName, ...
                        chanXName, dataH, dataX, timeH, timeX, transFnXtoH, segStartTime, ...
                        segEndTime, timeShift, samplFreq, logFid, debugLevel);
                case 'bilinear'

                    [crossCorrXHNew, crossCorrMaxXHNew, dataHClean] = bilinearcouplingcoeff(chanHName, ...
                        chanXName, dataH, dataX, timeH, timeX, segStartTime, segEndTime, ...
                        timeShift, samplFreq, logFid, 1, debugLevel);

                otherwise  

                    error('ERROR: Unknown coupling model.');
                end

				% construct an H data vector with the glitch replaced by the Hclean
                datSegIdxX = find(timeX >= segStartTime & timeX <= segEndTime);
                datSegIdxH = find(timeH >= segStartTime & timeH <= segEndTime);
				dataH2 = cell2mat(dataH);
				dataH2(datSegIdxH) = dataHClean;

				Data(coincIndex).dataH_clean = dataH2;


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % plot time-series data, if necessary
                if debugLevel >= 2

                    close all

                    trigCentTimeIdxH = find(timeH >= trigHCentTime, 1);
                    trigCentTimeIdxX = find(timeX >= trigXCentTime, 1);
                    t0X = min(timeX);
                    t0H = min(timeH);

                    if length(datSegIdxX) < 2 | length(datSegIdxH) < 2
                        error(['ERROR: length of the data segment is zero:',...
                        'len(X) =', num2str(length(datSegIdxX)), 'len(H) = ', ...
                        num2str(length(datSegIdxH))]);
                    end
                    
                    fh1 = plotdata(timeX-t0X, dataX{1}, 10, 'lin', 'lin', ...
                        ['t[sec] since ', num2str(t0X)], 'x(t)', ...
                        ['Time-series data:' , strrep(strcat(chanXName{1:end}),'_','\_')], ...
                        {'raw data','high-passed'},{'r--','b--'});
                    plot(timeX(datSegIdxX(1))-t0X, 0, 'k<', ...
                        timeX(datSegIdxX(end))-t0X, 0, 'k<')
                    plot(timeX(trigCentTimeIdxX)-t0X, 0, 'kx')
					xaxis(wreadStartTime-t0H+1, wreadEndTime-t0H-0.5);

                    fh2 = plotdata(timeH-t0H, dataH, 11, 'lin', 'lin', ...
                        ['t[sec] since ', num2str(t0H)], 'h(t)', ...
                        ['Time-series data:' , strrep(chanHName,'_','\_')], ...
                        {'raw data','high-passed'},{'c'});

                    fh2 = plotdata(timeH-t0H, dataH2, 11, 'lin', 'lin', ...
                        ['t[sec] since ', num2str(t0H)], 'h(t)', ...
                        ['Time-series data:' , strrep(chanHName,'_','\_'), ...
                        sprintf('\n r = %3.2f rNew = %3.2f', ...
                        crossCorrXHVec(coincIndex), crossCorrXHNew)], ...
                        {'h(t)','h_{clean}(t)'},{'b'}, [], 1);
                    plot(timeH(datSegIdxH(1))-t0H, 0, 'k<', ...
                        timeH(datSegIdxH(end))-t0H, 0, 'k<')
                    plot(timeH(trigCentTimeIdxH)-t0H, 0, 'kx')
					xaxis(wreadStartTime-t0H+1, wreadEndTime-t0H-0.5);

                    % save  plots
					wprintfig(fh1, [outDir '/plots/TimeSeries_', chanXName{1:end}, '_', ...
                        num2str(floor(t0X))]);
                    wprintfig(fh2, [outDir '/plots/TimeSeries_', chanHName, '_', ...
                        num2str(floor(t0H))]);
					
                end

            end     % END if statement checking whether the data read from framefiles is empty.
        end         % END if statement segStartTime > analysisStartTime ....
    end             % END for loop over coincident triggers 

	% save the mat file 
	save(outFile, 'Info', 'Data');
    
    % Close the file for writing the cross-correlation statistics for each
    % timeshift veto analysis.
else
    fprintf(logFid, ...
        'WARNING: No triggers with the specified correlation threshold...\n'); 
end

% close log file
fclose(logFid); 
