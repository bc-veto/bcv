function [rXH, rMaxXH, dataHClean] = bilinearcouplingcoeff(chanHName, ...
            chanXName, dataH, dataX, timeH, timeX, segStartTime, segEndTime, ...
            timeShift, samplFreq, logFid, subtractNoise, debugLevel)
% 
% BILINEARCOUPLINGCOEFF - calculate the cross-correlation coefficient 
% between the gravitational-wave channel H and the "projected" instrumental 
% channlel X. The noise in the instrumental channel X is "projected" 
% (transferred) to the domain of H using the data from another "slow" channel
% Y. 
% 
% i.e, projected channel X'(f) = X(f) Y(f)
% 
% 
% 
% Authors: P. Ajith, Aaron Pearlman
% 
% $Id: linearcouplingcoeff.m 242 2009-09-23 02:00:20Z ajith $
    
% Set the frequency range for the veto analysis.
MIN_FREQ = 40;
MAX_FREQ = 4000;

% Convert cell arrays to double vectors.
dataH = cell2mat(dataH);
dataX = cell2mat(dataX);     

% Check that the data segments from channels H and X are
% both non-empty.
if isempty(dataH) || isempty(dataX)

    fprintf(logFid, ...
        'ERROR: One or more data vectors are empty...\n');
    fprintf(logFid, ...
         'ERROR: length(dataH) = %d length(dataX) = %d\n', ...
         length(dataH), length(dataX));

% check that the data segments are equal 
elseif length(dataH) ~= length(dataX)

    fprintf(logFid, ...
        'ERROR: Different lengths. len(dataH) = %d len(dataX) = %d\n',...
            length(dataH), length(dataX));

else

    % Subtract the mean from the data in channels X and H.
    dataH = dataH - mean(dataH);
    dataX = dataX - mean(dataX);
    
    % Remove data that does not fall between segStartTime and
    % segEndTime in channel H.
    removeIndexChH = find(timeH < segStartTime | timeH >= segEndTime);
	%removeIndexChH = [];
    
    if ~isempty(removeIndexChH)
        timeH(removeIndexChH) = [];
        dataH(removeIndexChH) = [];
    end
                    
    % Remove data that does not fall between segStartTime and
    % segEndTime in channel X.
    removeIndexChX = find(timeX + timeShift < segStartTime ...
        | timeX + timeShift >= segEndTime);
	%removeIndexChX = [];
                    
    if ~isempty(removeIndexChX)
        timeX(removeIndexChX) = [];
        dataX(removeIndexChX) = [];
    end
                    
    % Set the parameters for calculating the FFT 
    nfft = length(dataX);
    %wind = hann(nfft);
    %wind = tukeywin(nfft);
    wind = ones(size(dataH));
                    
    % Calculate the FFT of the isolated data segment from channel H
    % and the product of X(t) and Y(t).

%    [fftChanH, freqVecH] = ...
%        spectrogram(dataH, wind, 0, nfft, samplFreq);
%    [fftXTransf, freqVecX] = ...
%        spectrogram(dataX, wind, 0, nfft, samplFreq);

    [fftChanH, freqVecH, ttH] = ...
        mSpecgram(dataH, wind, 0, nfft, samplFreq);
    [fftXTransf, freqVecX, ttX] = ...
        mSpecgram(dataX, wind, 0, nfft, samplFreq);

    % Check that the spectrogram channel H frequency vector
    % is the same size as the spectrogram channel X
    % frequency vector.
    if (size(freqVecH) ~= size(freqVecX)) 
        fprintf(logFid, ...
            'ERROR: Unequal size for freq. vectors.\n');
        fprintf(logFid, ...
            'ERROR: len(freqVecH) = %d len(freqVecX) = %d \n',...
                length(freqVecH), length(freqVecX));
        return;
    else
    
		% select the frequency band
        freqBandIdx = find(freqVecH > MIN_FREQ & freqVecH < MAX_FREQ);
                                
        % Calculate the cross-correlation statistic for the
        % segment of data in channels H and the "projected X" 
        [rXH, rMaxXH] = calccrosscorr(fftXTransf(freqBandIdx), fftChanH(freqBandIdx));

        % subtract the transferred noise from H, if necessary
        % also estimate the transfer function
        if subtractNoise == 1

           	transFnXYH =  fftXTransf.*conj(fftChanH)./(abs(fftXTransf).^2);

           	projXTransfH = fftXTransf.*transFnXYH;
           	projXTransfH = projectvector(projXTransfH, fftChanH);

           	%fftChanHClean = fftChanH-projXTransfH;
           	fftChanHClean = fftChanH;

			freqVecHNeg = flipud(-freqVecH(2:end-1));
			fftChanHCleanNeg = flipud(conj(fftChanHClean(2:end-1)));

			freqVecH2 = [freqVecH; freqVecHNeg];
			fftChanHClean = [fftChanHClean; fftChanHCleanNeg];

           	dataHClean = ifft(fftChanHClean);

        end

        % plot the spectra of channel H, X and X', if necessary
        if debugLevel >= 2
            fh1 = plotdata(freqVecX(freqBandIdx), abs(fftXTransf(freqBandIdx)), 20, 'log', 'log', 'f[Hz]', '|x(f)|', ...
                    ['Amp. spec:' , strrep(chanXName,'_','\_')], ...
                    {'x(t) y(t)'},{'r'});
            fh2 = plotdata(freqVecH(freqBandIdx), abs(fftChanH(freqBandIdx)), 21, 'log', 'log', 'f[Hz]', '|h(f)|', ...
                    ['Amp. spectra:' , strrep(chanHName,'_','\_')], ...
                    {'H','Xprime'},{'g'});
            fh2 = plotdata(freqVecH(freqBandIdx), abs(fftXTransf(freqBandIdx)), 21, 'log', 'log', 'f[Hz]', '|h(f)|', ...
                    strrep(['Amp. spec: H = ' , chanHName, ', X = ', chanXName{1}, ', Y = ', chanXName{2}],'_','\_'),...
                    {'H','Xprime'},{'b--'});
            if  subtractNoise == 1
                fh2 = plotdata(freqVecH2(freqBandIdx), abs(fftChanHClean(freqBandIdx)), 21, 'log', 'log', 'f[Hz]', '|h(f)|', ...
                    strrep(['Amp. spec: H = ' , chanHName, ', X = ', chanXName{1}, ', Y = ', chanXName{2}],'_','\_'),...
                    {'H','X_{prime}','H_{clean}'},{'k--'});

                fh3 = plotdata(freqVecH(freqBandIdx), abs(transFnXYH(freqBandIdx)), 22, 'log', 'log', 'f', 'Mag.', ...
                    strrep(['Transfer fn from (', chanXName{1}, ',', chanXName{2},') to ', chanHName],'_','\_'),...
                    '', {'r'}, 211);
                fh3 = plotdata(freqVecH(freqBandIdx), angle(transFnXYH(freqBandIdx)), 22, 'log', 'lin', 'f', 'Phase', ...
                    '', '', {'r'}, 212);
                saveas(fh3, ['plots/TransFn_', chanXName{1:end}, '_', chanHName, '_'...
			        num2str(floor(min(timeH)))], 'png');
            end

		    % save  plots
            saveas(fh1, ['plots/AmplSpectra_', chanXName{1:end}, '_', ...
			    num2str(floor(min(timeX)))], 'png');
		    saveas(fh2, ['plots/AmplSpectra_', chanHName, '_', ...
                num2str(floor(min(timeH)))], 'png');

        end

    end
end
                    

