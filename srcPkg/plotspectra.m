% DEBUGPLOTSPECTRA - Create amplitude spectra, power spectral densities,
% and time-frequency plots for debugging the vetoanalysis pipeline
% (channels H and X, along with the projected noise X').
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 06-08-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Vetoes for Binary-Coalescence Triggers Using Known Instrumental     %
%                                Couplings                                %
%                                                                         %
%                      Modified By: Aaron B. Pearlman                     %
%                        Mentor: Ajith Parameswaran                       %
%                        Creation Date: 06-08-2009                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
% 5. Plot the amplitude spectra of the channel H and projected X' data.

figure(3)
subplot(2, 1, 1)
loglog(freqChannelH, abs(specChannelH), 'r')
hold on
loglog(freqChannelX, abs(xPrime), 'k--')
grid on
title(sprintf('%s - Channel H and Projected X'' - Amplitude Spectra', ...
    chanHNameString))
xlabel('Frequency [Hz]')
ylabel('Abs[H(f)] OR Abs[X''(f)]')
xlim([40 1000])
legend('Channel H Spectra', 'Projected X'' Spectra')

%**************************************************************************
% 6. Plot the amplitude spectra of the channel X data.

figure(3)
subplot(2, 1, 2)
loglog(freqChannelX, abs(specChannelX), 'b')
grid on
title(sprintf('%s - Channel X - Amplitude Spectra', chanXNameString))
xlabel('Frequency [Hz]')
ylabel('Abs[X(f)]')
xlim([40 1000])

saveas(gcf, sprintf('%sAmplitudeSpectra.fig', plotString), 'fig')
saveas(gcf, sprintf('%sAmplitudeSpectra.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-AmplitudeSpectra.ps'

%**************************************************************************
% 7. Plot the ratio of the projected X' and channel H amplitude spectra,
% i.e. X' / H.

figure(4)
loglog(freqChannelH, abs(xPrime) / abs(specChannelH), 'r')
grid on
title(sprintf('%s - Channel X'' / H - Amplitude Spectra', chanXNameString))
xlabel('Frequency [Hz]')
ylabel('Abs[X''(f)] / Abs[H(f)]')
xlim([40 1000])

saveas(gcf, sprintf('%sSpectraRatioXPrimeToH.fig', plotString), 'fig')
saveas(gcf, sprintf('%sSpectraRatioXPrimeToH.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-SpectraRatioXPrimeToH.ps'

%**************************************************************************
% 8. Plot the magnitude of the transfer function before and after
% interpolation.

figure(5)
subplot(2, 1, 1)
loglog(transFnXtoH.frequency, abs(transFnXtoH.Txh), 'r-')
hold on
loglog(tfFreqIntp, tfMagIntp, 'b--')
grid on
title(sprintf('%s - Magnitude of Transfer Function (T_x_h)', ...
    chanXNameString))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
xlim([40 1000])
legend('T_x_h - Before', 'T_x_h - Interpolation')

%**************************************************************************
% 9. Plot the phase of the transfer function before and after
% interpolation.

figure(5)
subplot(2, 1, 2)
semilogx(transFnXtoH.frequency, unwrap(angle(transFnXtoH.Txh)), 'r-')
hold on
semilogx(tfFreqIntp, tfPhaseIntp, 'b--')
grid on
title(sprintf('%s - Phase of Transfer Function (T_x_h)', chanXNameString))
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
xlim([40 1000])
legend('T_x_h - Before', 'T_x_h - Interpolation')

saveas(gcf, sprintf('%sTransferFunction.fig', plotString), 'fig')
saveas(gcf, sprintf('%sTransferFunction.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-TransferFunction.ps'

%**************************************************************************
% 10. Calculate the mean of the channel H and X' amplitude spectra vectors
% over all frequencies.
sampleValue = mean(abs(specChannelH) ./ abs(xPrime));

%**************************************************************************
% 11. Compute and plot the average specta for channel H before and after
% high-passing.

nfft = length(rawDataChH);
nOverlap = nfft/2;
wind = hann(nfft);
fs = rawSampleFreq;

[PxxHRaw, FRawChH] = pwelch(rawDataChH, wind, nOverlap, nfft, fs);
[PxxHHP, FHighPassChH] = pwelch(rawDataChHHP, wind, nOverlap, nfft, fs);

% Remove frequencies from the power spectral densities of channel H where
% the transfer function is not defined.
spectraDataRemoveIndexChH = find(FRawChH < TxhFreqMin ...
    | FRawChH > TxhFreqMax);

if length(FRawChH) == length(FHighPassChH) ...
        && length(PxxHRaw) == length(PxxHHP) ...
        && ~isempty(spectraDataRemoveIndexChH)
    
    FRawChH(spectraDataRemoveIndexChH) = [];
    PxxHRaw(spectraDataRemoveIndexChH) = [];
    
    FHighPassChH(spectraDataRemoveIndexChH) = [];
    PxxHHP(spectraDataRemoveIndexChH) = [];
end

figure(6)
loglog(FRawChH, sqrt(PxxHRaw), 'r')
hold on
loglog(FHighPassChH, sqrt(PxxHHP), 'c')

% Read in 64 seconds of time-series data from channel H containing no
% glitches to use as a background.
backgroundStartTimeChH = 928888265;
backgroundEndTimeChH = 928888329;

[backgroundDataChH, backgroundSampleFrequenciesChH] = ...
    wreaddata(frameCache, chanHName, frameTypeChanH, ...
    backgroundStartTimeChH, backgroundEndTimeChH, 0, debugLevel);

% Subtract the mean from the background data in channel H.
backgroundDataChH = cell2mat(backgroundDataChH);
backgroundDataChH = backgroundDataChH - mean(backgroundDataChH);

% High-pass the background data in channel H above 40 Hz.
tiling.backgroundSampleFrequenciesChH = backgroundSampleFrequenciesChH(1);
tiling.backgroundHighPassCutoffChH = 40;

[highPassedBackgroundDataChH, coefficientsChH] = ...
    whighpass(backgroundDataChH, tiling);

% Subtract the mean from the high-passed background data in channel H.
highPassedBackgroundDataChH = cell2mat(highPassedBackgroundDataChH);
highPassedBackgroundDataChH = ...
    highPassedBackgroundDataChH - mean(highPassedBackgroundDataChH);

% Compute the average spectra of the background data in channel H.
nfft = 8192;
nOverlap = nfft/2;
wind = hann(nfft);
fs = backgroundSampleFrequenciesChH;

[PxxHBackRaw, FHBackRaw] = pwelch(backgroundDataChH, wind, nOverlap, ...
    nfft, fs);
[PxxHBackHP, FHBackHP] = pwelch(highPassedBackgroundDataChH, wind, ...
    nOverlap, nfft, fs);

% Remove frequencies from the background power spectral densities of
% channel H where the transfer function is not defined.
spectraBackDataRemoveIndexChH = find(FHBackRaw < TxhFreqMin ...
    | FHBackRaw > TxhFreqMax);

if length(FHBackRaw) == length(FHBackHP) ...
        && length(PxxHBackRaw) == length(PxxHBackHP) ...
        && ~isempty(spectraBackDataRemoveIndexChH)
    
    FHBackRaw(spectraBackDataRemoveIndexChH) = [];
    PxxHBackRaw(spectraBackDataRemoveIndexChH) = [];
    
    FHBackHP(spectraBackDataRemoveIndexChH) = [];
    PxxHBackHP(spectraBackDataRemoveIndexChH) = [];
end

hold on
loglog(FHBackRaw, sqrt(PxxHBackRaw), 'b')
hold on
loglog(FHBackHP, sqrt(PxxHBackHP), 'k')

%**************************************************************************
% 12. Compute and plot the average specta for channel X before and after
% high-passing.

nfft = length(rawDataChX);
nOverlap = nfft/2;
wind = hann(nfft);
fs = rawSampleFreq;

[PxxXRaw, FRawChX] = pwelch(rawDataChX, wind, nOverlap, nfft, fs);
[PxxXHP, FHighPassChX] = pwelch(rawDataChXHP, wind, nOverlap, nfft, fs);

% Remove frequencies from the power spectral densities of channel X where
% the transfer function is not defined.
spectraDataRemoveIndexChX = find(FRawChX < TxhFreqMin ...
    | FRawChX > TxhFreqMax);

if length(FRawChX) == length(FHighPassChX) ...
        && length(PxxXRaw) == length(PxxXHP) ...
        && ~isempty(spectraDataRemoveIndexChX)
    
    FRawChX(spectraDataRemoveIndexChX) = [];
    PxxXRaw(spectraDataRemoveIndexChX) = [];
    
    FHighPassChX(spectraDataRemoveIndexChX) = [];
    PxxXHP(spectraDataRemoveIndexChX) = [];
end

figure(7)
loglog(FRawChX, sqrt(PxxXRaw), 'r')
hold on
loglog(FHighPassChX, sqrt(PxxXHP), 'c')

% Read in 64 seconds of time-series data from channel X containing no
% glitches to use as a background.
backgroundStartTimeChX = 928888800;
backgroundEndTimeChX = 928888864;

[backgroundDataChX, backgroundSampleFrequenciesChX] = ...
    wreaddata(frameCache, chanXName, frameTypeChanX, ...
    backgroundStartTimeChX, backgroundEndTimeChX, 0, debugLevel);

% Subtract the mean from the background data in channel X.
backgroundDataChX = cell2mat(backgroundDataChX);
backgroundDataChX = backgroundDataChX - mean(backgroundDataChX);

% High-pass the background data in channel X above 40 Hz.
tiling.backgroundSampleFrequenciesChX = backgroundSampleFrequenciesChX(1);
tiling.backgroundHighPassCutoffChX = 40;

[highPassedBackgroundDataChX, coefficientsChX] = ...
    whighpass(backgroundDataChX, tiling);

% Subtract the mean from the high-passed background data in channel X.
highPassedBackgroundDataChX = cell2mat(highPassedBackgroundDataChX);
highPassedBackgroundDataChX = ...
    highPassedBackgroundDataChX - mean(highPassedBackgroundDataChX);

% Compute the average spectra of the background data in channel X.
nfft = 8192;
nOverlap = nfft/2;
wind = hann(nfft);
fs = backgroundSampleFrequenciesChX;

[PxxXBackRaw, FXBackRaw] = pwelch(backgroundDataChX, wind, nOverlap, ...
    nfft, fs);
[PxxXBackHP, FXBackHP] = pwelch(highPassedBackgroundDataChX, wind, ...
    nOverlap, nfft, fs);

% Remove frequencies from the background power spectral densities of
% channel X where the transfer function is not defined.
spectraBackDataRemoveIndexChX = find(FXBackRaw < TxhFreqMin ...
    | FXBackRaw > TxhFreqMax);

if length(FXBackRaw) == length(FXBackHP) ...
        && length(PxxXBackRaw) == length(PxxXBackHP) ...
        && ~isempty(spectraBackDataRemoveIndexChX)
    
    FXBackRaw(spectraBackDataRemoveIndexChX) = [];
    PxxXBackRaw(spectraBackDataRemoveIndexChX) = [];
    
    FXBackHP(spectraBackDataRemoveIndexChX) = [];
    PxxXBackHP(spectraBackDataRemoveIndexChX) = [];
end

hold on
loglog(FXBackRaw, sqrt(PxxXBackRaw), 'b')
hold on
loglog(FXBackHP, sqrt(PxxXBackHP), 'k')
grid on
title(sprintf('%s - Channel X - Power Spectral Densities of the Raw and High-Passed Data', ...
    chanXNameString))
xlabel('Frequency [Hz]')
ylabel('X(f) / Sqrt[Hz]')
xlim([40 1000])
legend('Raw Data Channel X', 'High-Passed Data Channel X', ...
    'Background - Raw Data Channel X', ...
    'Background - High-Passed Data Channel X')

saveas(gcf, sprintf('%sPowerSpectralDensityChX.fig', plotString), 'fig')
saveas(gcf, sprintf('%sPowerSpectralDensityChX.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-PowerSpectralDensityChX.ps'

%**************************************************************************
% 13. Compute and plot the average spectra for the projected noise X'.

% Calculate the frequency resolution for interpolating the transfer
% function from channel X to H.
frequencyResolutionTxh = FHighPassChX(2) - FHighPassChX(1);

% Interpolate the transfer function with the frequency resolution from
% channel X.
[tfFreqIntp, tfMagIntp, tfPhaseIntp] = ...
    interpolatetransfn(transFnXtoH.frequency, abs(transFnXtoH.Txh), ...
    unwrap(angle(transFnXtoH.Txh)), frequencyResolutionTxh);

TxhInterp = tfMagIntp .* exp(1i * tfPhaseIntp);

% X(~)' = X(~) * Txh.Txh (interpolated transfer function).
% The transfer function is in units of displacement. We have to divide this
% by the length of the interferometer arms to convert into units of strain.
% For H1 and L1, the length is 4000 m, and for H2, the length is 2000 m.
xPrimeSpectra = sqrt(PxxXHP) .* abs(TxhInterp / 4000);

% Compute and plot the X' specta.
figure(6)
hold on
loglog(FHighPassChX, xPrimeSpectra, 'm')
grid on
title(sprintf('%s - Channel H - Power Spectral Densities of the Raw and High-Passed Data', ...
    chanHNameString))
xlabel('Frequency [Hz]')
ylabel('H(f) / Sqrt[Hz]')
xlim([40 1000])
legend('Raw Data Channel H', 'High-Passed Data Channel H', ...
    'Background - Raw Data Channel H', ...
    'Background - High-Passed Data Channel H', 'Projected Noise X''')

saveas(gcf, sprintf('%sPowerSpectralDensityChH.fig', plotString), 'fig')
saveas(gcf, sprintf('%sPowerSpectralDensityChH.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-PowerSpectralDensityChH.ps'

%**************************************************************************
% 14. Create a time-frequency plot of the raw data from channel H.

% Compute the specta before and after high-passing the data.
if length(rawDataChH) == length(rawDataChHHP)
    nfft = length(rawDataChH);
end

nOverlap = nfft/2;
wind = hann(nfft);
fs = rawSampleFreq(1);

[specChannelHRaw, FChHRaw, TChHRaw] = spectrogram(rawDataChH, wind, ...
    nOverlap, nfft, fs);
[specChannelHHP, FChHHP, TChHP] = spectrogram(rawDataChHHP, wind, ...
    nOverlap, nfft, fs);

% Create a time-frequency plot of the raw data.
figure(8)
imagesc(TChHRaw, FChHRaw, log10(abs(specChannelHRaw)))
grid on
colorbar
title(sprintf('Time-Frequency Plot of the Raw Data - Channel H: %s', ...
    chanHNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Frequency [Hz]')
ylim([40 1000])

saveas(gcf, sprintf('%sRawDataTimeFreqChH.fig', plotString), 'fig')
saveas(gcf, sprintf('%sRawDataTimeFreqChH.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-RawDataTimeFreqChH.ps'

%**************************************************************************
% 15. Create a time-frequency plot of the high-passed data from channel H.

figure(9)
imagesc(TChHP, FChHHP, log10(abs(specChannelHHP)))
grid on
colorbar
title(sprintf('Time-Frequency Plot of the High-Passed Data - Channel H: %s', ...
    chanHNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Frequency [Hz]')
ylim([40 1000])

saveas(gcf, sprintf('%sHighPassedTimeFreqChH.fig', plotString), 'fig')
saveas(gcf, sprintf('%sHighPassedTimeFreqChH.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-HighPassedTimeFreqChH.ps'

%**************************************************************************
% 16. Create a time-frequency plot of the raw data from channel X.

% Compute the specta before and after high-passing the data.
if length(rawDataChX) == length(rawDataChXHP)
    nfft = length(rawDataChX);
end

nOverlap = nfft/2;
wind = hann(nfft);
fs = rawSampleFreq(1);

[specChannelXRaw, FChXRaw, TChXRaw] = spectrogram(rawDataChX, wind, ...
    nOverlap, nfft, fs);
[specChannelXHP, FChXHP, TChXP] = spectrogram(rawDataChXHP, wind, ...
    nOverlap, nfft, fs);

% Create a time-frequency plot of the raw data.
figure(10)
imagesc(TChXRaw, FChXRaw, log10(abs(specChannelXRaw)))
grid on
colorbar
title(sprintf('Time-Frequency Plot of the Raw Data - Channel X: %s', ...
    chanXNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Frequency [Hz]')
ylim([40 1000])

saveas(gcf, sprintf('%sRawDataTimeFreqChX.fig', plotString), 'fig')
saveas(gcf, sprintf('%sRawDataTimeFreqChX.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-RawDataTimeFreqChX.ps'

%**************************************************************************
% 17. Create a time-frequency plot of the high-passed data from channel X.

figure(11)
imagesc(TChXP, FChXHP, log10(abs(specChannelXHP)))
grid on
colorbar
title(sprintf('Time-Frequency Plot of the High-Passed Data - Channel X: %s', ...
    chanXNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Frequency [Hz]')
ylim([40 1000])

saveas(gcf, sprintf('%sHighPassedTimeFreqChX.fig', plotString), 'fig')
saveas(gcf, sprintf('%sHighPassedTimeFreqChX.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-HighPassedTimeFreqChH.ps'
return
%**************************************************************************