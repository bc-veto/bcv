clear 

tic 

addpath('/archive/home/ajith/working/dc/veto/omega_veto/src');

dbFile = '/archive/home/ajith/working/dc/veto/omega_veto/test/noisesub_config/H1-LSC-MICH_CTRL+H1-ASC-QPDX_P_H1-LSC-MICH_CTRL+H1-ASC-QPDX_P-data.db';
frameCacheFile = '/archive/home/ajith/working/dc/veto/omega_veto/test/noisesub_config/LHO.cache';
chanHName = 'H1:LDAS-STRAIN';
chanXName = {'H1:LSC-MICH_CTRL','H1:ASC-QPDX_P'};
frameTypeChanH = 'H1_LDAS_C02_L2';
frameTypeChanX = {'H1_RDS_R_L1','H1_RDS_R_L1'};
samplFreqH = 4096;
samplFreqX = 4096;
highPassCutoff = 32;
couplingModel = 'bilinear';
transFnXtoH = 'null';
analysisStartTime = 957917560.0;
analysisEndTime = 959114559.0;
corrThreshMin = 0.7;
corrThreshMax = 1.0;
outDir = '2010-07-17-NoiseSubTest';
logFile = 'test.log';
outFile = 'Test.mat';
debugLevel = 1;

subtractnoise(dbFile, frameCacheFile, chanHName, chanXName, frameTypeChanH, ...
 	frameTypeChanX, samplFreqH, samplFreqX, highPassCutoff,  ...
 	couplingModel, transFnXtoH, analysisStartTime, analysisEndTime, ...
 	corrThreshMin, corrThreshMax, outDir, logFile, outFile, debugLevel)

toc

