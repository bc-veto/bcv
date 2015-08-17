#!/usr/bin/python
import os;
import sys;
import math;
import string;

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

cmd = './omegaveto %s %s %s %s %2.1f %2.1f %2.1f %s %s %2.1f %2.1f %d %d' %(segmentFile, \
		configurationFile, frameCacheFile, couplingModel, highPassCutoff, \
		trigSignThreshX, trigSignThreshH, outDir, logFile, timeShiftMin, \
		timeShiftMax, numTimeShifts, debugLevel);

os.system('source MatlabSetup_R2008a_glnxa64.sh');
os.system(cmd);

