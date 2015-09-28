function [rXH, rMaxXH] = linearcouplingcoeff(chanHName, chanXName, ...
            dataH, dataX, timeH, timeX, transFnXtoH, segStartTime, ...
            segEndTime, timeShift, samplFreq, logFid, debugLevel)
% 
% LINEARCOUPLINGCOEFF - calculate the cross-correlation coefficient 
% between the gravitational-wave channel H and the "projected" instrumental 
% channlel X. The noise in the instrumental channel X is "projected" 
% (transferred) to the domain of H using a linear-coupling transfer 
% function Txh. 
% 
% i.e, projected channel X'(f) = X(f) Txh(f)
% 
% 
% 
% Authors: Aaron Pearlman, P. Ajith 
% 
% $Id: linearcouplingcoeff.m 277 2009-10-31 00:02:03Z ajith $

% Set the frequency range for the veto analysis.
MIN_FREQ = 10;      % min freq of the analysis band (Hz)
MAX_FREQ = 4000;    % max freq of the analysis band (Hz)
IFO_LENGTH = 4000;  % length of the interferometr (meters)

% Convert cell arrays to double vectors.
dataH = cell2mat(dataH);
dataX = cell2mat(dataX);

% Check that the data segments from channels H and X are
% both non-empty.
if isempty(dataH) || isempty(dataX)

    fprintf(logFid, ...
        'ERROR: One or more data vectors are empty...\n');
    fprintf(logFid, ...
         'ERROR: length(dataH) = %d length(dataX) = %d...\n', ...
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
    removeIndexChH = find(timeH < segStartTime ...
        | timeH >= segEndTime);
    
    if ~isempty(removeIndexChH)
        timeH(removeIndexChH) = [];
        dataH(removeIndexChH) = [];
    end
                    
    % Remove data that does not fall between segStartTime and
    % segEndTime in channel X.
    removeIndexChX = find(timeX + timeShift < segStartTime ...
        | timeX + timeShift >= segEndTime);
                    
    if ~isempty(removeIndexChX)
        timeX(removeIndexChX) = [];
        dataX(removeIndexChX) = [];
    end
                    
    % Set the parameters for calculating the FFT of the isolated data 
    % segment from channels H and X.
    nfft = length(dataX);
    wind = hann(nfft);
           
    % Calculate the FFT of the isolated data segment from
    % channel H.
    [fftChanH, freqVecH] = ...
        spectrogram(dataH, wind, 0, nfft, samplFreq);
                        
    % Calculate the amplitude spectra of the isolated data
    % segment from channel X.
    [fftChanX, freqVecX] = ...
        spectrogram(dataX, wind, 0, nfft, samplFreq);
                        
    % Remove data from the transfer function that does not
    % lie within the desired frequency range.
    specFreqRemovIdx = find(transFnXtoH.frequency < MIN_FREQ | ...
        transFnXtoH.frequency > MAX_FREQ);
                        
    if ~isempty(specFreqRemovIdx)
        transFnXtoH.frequency(specFreqRemovIdx) = [];
        transFnXtoH.Txh(specFreqRemovIdx) = [];
    end
                        
    % Find the minimum and maximum frequency of the
    % transfer function.
    TxhFreqMin = min(transFnXtoH.frequency);
    TxhFreqMax = max(transFnXtoH.frequency);
                        
    % Check that the spectrogram channel H frequency vector
    % is the same size as the spectrogram channel X
    % frequency vector.
    if (size(freqVecH) ~= size(freqVecX))
        fprintf(logFid, ...
            'ERROR: Unequal size for freq. vectors.\n');
        fprintf(logFid, ...
            'ERROR: len(freqVecH) = %d len(freqVecX) = %d.\n',...
                length(freqVecH), length(freqVecX));
        return;
    else
    
        % Remove spectrogram data from channels H and X
        % that do not lie within the desired frequency
        % range of the transfer function.
        frequencySpectraRemoveIndex = ...
            find(freqVecH < TxhFreqMin ...
            | freqVecH > TxhFreqMax);
                                
        if ~isempty(frequencySpectraRemoveIndex)
            freqVecH(frequencySpectraRemoveIndex) = [];
            fftChanH(frequencySpectraRemoveIndex) = [];
            freqVecX(frequencySpectraRemoveIndex) = [];
            fftChanX(frequencySpectraRemoveIndex) = [];
        end
                                
        % Calculate the frequency resolution for
        % interpolating the transfer function from channel
        % X to H.
        freqResolTxh = samplFreq / length(dataX);
                                
        % Interpolate the transfer function with the
        % frequency resolution from channel X.
        [tfFreqIntp, tfMagIntp, tfPhaseIntp] = interpolatetransfn(...
            transFnXtoH.frequency, abs(transFnXtoH.Txh), ...
            unwrap(angle(transFnXtoH.Txh)), freqResolTxh);
                                
        TxhInterp = tfMagIntp .* exp(1i * tfPhaseIntp);
                                
        % X(~)' = X(~) * Txh (interpolated transfer
        % function).
        %
        % The transfer function is in units of
        % displacement.
        %
        % We have to divide this by the length of the
        % interferometer arms to convert into units of
        % strain.
        %
        % For H1 and L1, the length is 4000 m, and for H2,
        % the length is 2000 m.
        if size(fftChanX) == size(TxhInterp)
            xPrime = fftChanX .* (TxhInterp / IFO_LENGTH);
        else
            fprintf(logFid, ...
                'ERROR: size(fftChanX) = %d size(TxhInterp) = %d\n', ...
                size(fftChanX), size(TxhInterp));
            fprintf(logFid, ...
                'ERROR: fftChanX and TxhInterp have different sizes.\n');
        end
                                
        % Calculate the cross-correlation statistic for the
        % segment of data in channels H and X' (where the "projected
        % data" X' = X Txh)
        [rXH, rMaxXH] = calccrosscorr(xPrime, fftChanH);

        % plot the spectra of channel H, X and X', if necessary
        if debugLevel >= 2
            fh1 = plotdata(freqVecX, abs(fftChanX), 20, 'log', 'log', 'f[Hz]', '|x(f)|', ...
                    ['Ampl. spectra:' , strrep(chanXName,'_','\_')], ...
                    {'X'},{'r'});
            fh2 = plotdata(freqVecH, abs(fftChanH), 21, 'log', 'log', 'f[Hz]', '|h(f)|', ...
                    ['Ampl. spectra:' , strrep(chanHName,'_','\_')], ...
                    {'H','Projected X'},{'b'});
            fh2 = plotdata(freqVecH, abs(xPrime), 21, 'log', 'log', 'f[Hz]', '|h(f)|', ...
                    strrep(['Ampl. spectra: H = ' , chanHName, ', X = ', chanXName],'_','\_'), ...
                    {'H','Projected X'},{'k--'});

            % save  plots
            saveas(fh1, ['AmpSpectrum', chanXName{1:end}, '_', ...
                num2str(floor(min(timeX)))], 'png')
            saveas(fh2, ['AmpSpectrum', chanHName, '_', ...
                num2str(floor(min(timeH)))], 'png')

        end

    end
end
                    

