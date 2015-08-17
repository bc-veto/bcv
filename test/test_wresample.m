
clear 
addpath('../src');

cacheFile = '/archive/home/ajith/working/dc/veto/omega-veto/framecache/E14_nodes.txt';
frameCache = loadframecache(cacheFile);
startTime = utc2gps('2009-06-12 21:57:00');
stopTime = utc2gps('2009-06-12 21:57:34');
debugLevel = 2;

channelName = 'H1:LDAS-STRAIN';
frameTypes = 'H1_LDAS_C00_L2';

%channelName = 'H1:LSC-MICH_CTRL';
%frameTypes = 'R';

%channelName = 'H1:ASC-QPDX_P';
%frameTypes = 'R';

samplFreq2 = 4096;

% Read the data.
[data, samplFreq] = ...
      wreaddata(frameCache, channelName, frameTypes, startTime, ...
      stopTime, 0, debugLevel);

samplFreq

% high pass and resample the data 
tiling.sampleFrequency = samplFreq;
tiling.highPassCutoff = 32;
[data_hp] = whighpass(data, tiling);
data_hp_res = wresample(data_hp, samplFreq, samplFreq2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot time domain data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time  = [startTime+1: 1/samplFreq  : (stopTime-1) - 1/samplFreq]  - startTime;
time2 = [startTime+1: 1/samplFreq2 : (stopTime-1) - 1/samplFreq2] - startTime;

data = cell2mat(data);
data_hp = cell2mat(data_hp);
data_hp_res = cell2mat(data_hp_res);

% throw away the first and last second - junk due to high-pass and resampling 
data = data(samplFreq+1:end-samplFreq);
data_hp = data_hp(samplFreq+1:end-samplFreq);
data_hp_res = data_hp_res(samplFreq2+1:end-samplFreq2);

data = data - mean(data);
data_hp = data_hp - mean(data_hp);
data_hp_res = data_hp_res - mean(data_hp_res);

fh1 = figure(1)
subplot(211)
plot(time, data, 'r')
ylabel('Data')
hold on
grid on 
title(sprintf('%s - Channel X - Time Series', channelName))
subplot(212)
plot(time, data_hp, 'r', time2, data_hp_res, 'k--')
grid on
xlabel(['Time Since ', num2str(startTime) ' [s]'])
ylabel('Data')
legend('High passed','Resampled')
saveas(fh1, ['TimeSeries_' channelName '.fig'], 'fig')
saveas(fh1, ['TimeSeries_' channelName '.pdf'], 'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot FFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft = length(data)
nfft_res = length(data_hp_res)
wind = ones(size(data));
wind_res = ones(size(data_hp_res));

[fftData, f, t] = mSpecgram(data, wind, 0, nfft, samplFreq);
[fftData_hp, f, t] = mSpecgram(data_hp, wind, 0, nfft, samplFreq);
[fftData_hp_res, f2, t] = mSpecgram(data_hp_res, wind_res, 0, nfft_res, samplFreq2);

fh2 = figure(2)
subplot(121)
semilogx(f,real(fftData),'r', f, real(fftData_hp),'k--', f2, real(fftData_hp_res),'g--')
grid on
xlabel('f')
ylabel('Real part')
subplot(122)
semilogx(f,imag(fftData),'r', f, imag(fftData_hp),'k--', f2, imag(fftData_hp_res),'g--')
grid on
xlabel('f')
ylabel('Im part')
legend('raw','high-passed','resampled')
allxaxis(40, 2200)

saveas(fh2, ['FFT_ReIm_' channelName '.fig'], 'fig')
saveas(fh2, ['FFT_ReIm_' channelName '.pdf'], 'pdf')

fh3 = figure(3)
subplot(121)
loglog(f,abs(fftData)/nfft,'r', f, abs(fftData_hp)/nfft,'k--', f2, abs(fftData_hp_res)/nfft_res, 'g--') 
xlabel('f')
ylabel('Magnitude')
grid on
subplot(122)
semilogx(f,(angle(fftData)),'r', f, (angle(fftData_hp)),'k--', f2, (angle(fftData_hp_res)),'g--')
xlabel('f')
grid on
ylabel('Phase')
legend('raw','high-passed','resampled')

saveas(fh3, ['FFT_MagPhase_' channelName '.fig'], 'fig')
saveas(fh3, ['FFT_MagPhase_' channelName '.pdf'], 'pdf')
allxaxis(40, 2200)




