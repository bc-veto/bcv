clear 

tic 

addpath('/archive/home/ajith/Runs/BCV_S6/omega_veto/src');

analysisStartTime = utc2gps('2009-06-12 20:57:54');
analysisEndTime = utc2gps('2009-06-13 17:00:00');
generateReport = true;
debugLevel = 0;
highPassCutoff = 40;
trigSignThreshX = 100;
trigSignThreshH = 100;
couplingModel = 'bilinear';

timeShiftMin = 0;
timeShiftMax = 0;
numTimeShifts = 1;

configurationFile = '/archive/home/ajith/Runs/BCV_S6/omega_veto/test/test.conf';
frameCacheFile = '/archive/home/ajith/Runs/BCV_S6/omega_veto/test/E14.cache';
segmentFile = '/archive/home/ajith/Runs/BCV_S6/omega_veto/test/E14.seg';
outDir = 'TestResults';
logFile = 'test.log';

omegaveto(segmentFile, configurationFile, frameCacheFile, ...
        couplingModel, highPassCutoff, trigSignThreshX, trigSignThreshH, ...
        outDir, logFile, timeShiftMin, timeShiftMax, numTimeShifts, debugLevel);

toc

