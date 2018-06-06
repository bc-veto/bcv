#!/bin/env python
#/home/detchar/opt/gwpysoft/bin/python

'''
- bcv_runscript - perform veto analysis using different coupling models between 
- the instrumental channel and the gravitational-wave channel. 
- 
- usage: bcv_runscript.py segmentFile, configurationFile, frameCacheFile, ...
-        couplingModel, highPassCutoff, trigSignThreshX, trigSignThreshH, ...
-        outDir, logFile, timeShiftMin, timeShiftMax, numTimeShifts, debugLevel
-        
- configurationFile 	- Path name of channel configuration file
- frameCacheFileH 	- Path name of h(t) framecache file (Now lalcache)
- frameCacheFileX       - Path name of the x(t) and y(t) framecache file (Now lalcache)
- outDir 		- Directory to write results
- logFile		- string specifying the name of the log file 
- debugLevel 		- Verboseness of debug level output
-
- The configuration file is an ASCII text file describing the parameters
- for each channel to be analyzed. The entries should have the following
- syntax:
-
- We will assume that the first channel that is read from the configuration
- file is channel H (the h(t) detector output data vector), followed by a
- X channel and Y channel (the "fast" and "slow" aux. channles). 
-
- Aaron B. Pearlman <aaronp1@umbc.edu>, 22-06-09
- P. Ajith <ajith@caltech.edu>, 23-07-09
- Tomoki Isogai <isogait@carleton.edu>

- Re-written by Sudarshan Ghonge <sudu.ghonge@gmail.com> 27-09-15
 '''
import sys
import string
import numpy as np
import os
from time import strftime
import trigstruct
import vetoanalysis
import bcv
import gwpy.time as gtime
import re

# Class to store the configuration information of a particular channel
class ConfigurationClass:
  def __init__(self):
    self.channelName = []
    self.frameType = []
    self.sampleFrequency = []
    self.triggerListChH = []
    self.triggerListChX = []
    self.transferFunctionXtoH = []
# Class to store the transfer function of a particular fast channel
class TransferFunctionXtoH:
  def __init__(self, Freq, txh):
    self.frequency = Freq
    self.Txh = txh

# Function to check if the file exists AND is non-empty
def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

#if(len(sys.argv)<13):
  #print 'Veto analysis pipeline using instrumental coupling models\n'
  #print 'Usage (Python executable from unix shell)\n'
  #print '\tbcv_runscript.py <segmentFile> <configurationFile> <frameCacheFile>'
  #print '\t\t<couplingModel> <highPassCutoff> <trigSignThreshX> <trigSignThreshH>'
  #print '\t\t <outDir> <logFile> <timeShiftMin> <timeShiftMax> <numTimeShifts>'
  #print 'debugLevel'
  
  #print '\nSudarshan Ghonge, P Ajith 2015'
  #sys.exit()

import argparse
parser = argparse.ArgumentParser(description='Veto analysis pipeline using instrumental coupling models')
parser.add_argument('segmentFile', type=str, help='Name of file where information of the segments (in seconds) to be analysed is stored')
parser.add_argument('configurationFile', type=str, help='Name of the file which describes the configuration of the analysis to be done')
parser.add_argument('frameCacheFileH', type=str, help='Name of the file where information about h(t) frame files is stored')
parser.add_argument('frameCacheFileX', type=str, help='Name of the file where information about the x(t) and y(t) files is stored')
parser.add_argument('couplingModel', type=str, help='Choose between linear and bilinear', choices=['linear', 'bilinear'])
parser.add_argument('highPassCutoff', type=float, help='Provide high pass frequency cut-off')
parser.add_argument('trigSignThreshX', type=float, help='Provide minimum SNR value needed in X triggers')
parser.add_argument('trigSignThreshH', type=float, help='Provide minimum SNR value needed in H triggers')
parser.add_argument('outDir', type=str, help='Output directory of output files')
parser.add_argument('logFileName', type=str, help='Name of log file')
parser.add_argument('timeShiftMin', type=float, help='Minimum amount of time shift between H and X data')
parser.add_argument('timeShiftMax', type=float, help='Maximum amount of time shift between H and X data')
parser.add_argument('numTimeShifts', type=int, help='Number of time shifts between H and X data')
parser.add_argument('debugLevel', type=int,help='Provide debug level', choices=[0, 1, 2])
parser.add_argument('sampleFrequency', type=float, help='Sample Frequency')


args = parser.parse_args()

segmentFile = args.segmentFile
configurationFile = args.configurationFile
frameCacheFileH = args.frameCacheFileH
frameCacheFileX = args.frameCacheFileX
couplingModel = args.couplingModel
highPassCutoff = args.highPassCutoff
trigSignThreshX = args.trigSignThreshX
trigSignThreshH = args.trigSignThreshH
outDir = args.outDir
logFileName = args.logFileName
timeShiftMin = args.timeShiftMin
timeShiftMax = args.timeShiftMax
numTimeShifts = args.numTimeShifts
debugLevel = args.debugLevel
sampleFrequency = args.sampleFrequency


print "Reading segment file %s..\n" %(segmentFile)

# Read segment file
segList = np.loadtxt(segmentFile).reshape(-1, 4)

segStartTimeVec = segList[:, 1]
segEndTimeVec = segList[:, 2]
  
# Read configuration File

print 'Reading configuration file %s.....\n' %(configurationFile)

configuration = []
channelNumber = -1
sectionIndex = []
sectionName = []
sectionNumber = 0
sectionStart = []

# Open configuration file for reading
with open(configurationFile, 'r') as configurationFID:
  for line in configurationFID:
    #Remove comments
    commentIndex = line.find('#')
    if(commentIndex>=0):
      line = line[0:commentIndex]
      
    line = line.strip()
    
    if(len(line)==0):
      continue
    
    if(line[0]=='['):
      commaIndex = line.find(',')
    
      if(commaIndex<0):
        print 'Syntax error processing configuration file %s:\n%s\n'%(configurationFile,line)
      
      # Increment section number
      sectionNumber+=1
      #Extract section index
      msecIndex = line[1:commaIndex-1]
      #Extract section name
      msecName = line[commaIndex+1:len(line)-2]
      
      #Remove leading and trailing blanks
      msecIndex = msecIndex.strip()
      
      msecName = msecName.strip()
      
      #Standardize names
      msecIndex = msecIndex.replace(';', ':')
      msecName = msecName.replace(';', ':')
      
      msecStart = channelNumber+1
      
      sectionIndex.append(msecIndex)
      sectionName.append(msecName)
      sectionStart.append(msecStart)
      continue
    
    # Check for beginning of new channel config
    if(line=='{'):
      # Increment channel number
      channelNumber+=1
      # Allocate space for the config struct to store data for the
      # first two channels read from the config file.
      if(channelNumber<=2):
	# Initialize config parameters
	mconf = ConfigurationClass()
	configuration.append(mconf)
	
      continue
    # Check for end of existing channel configuration.
    if(line=='}'):
      # Validate channel config
      if(len(configuration[channelNumber].channelName)==0):
	sys.exit('Channel name not specified for channel number %d' %(channelNumber))

      if(len(configuration[channelNumber].frameType)==0):
	sys.exit('Frame not specified for channel number %d' %(channelNumber))
	
      if(configuration[channelNumber].sampleFrequency==None):
	print 'Warning: Incomplete configuration.\n Sample frequency not specified for channel number %d'%(channelNumber)
	
	configuration[channelNumber].sampleFrequency=sampleFrequency
      
      if((channelNumber==0)&(len(configuration[channelNumber].triggerListChH)==0)):
	print 'Warning: Incomplete configuration\nChannel H trigger list not specified for channel number %d'%(channelNumber)
	
	configuration[channelNumber].triggerListChH = defaultTriggerListChH

      if((channelNumber==1)&(len(configuration[channelNumber].triggerListChX)==0)):
	print channelNumber==1
	print 'Warning: Incomplete configuration\nChannel X trigger list not specified for channel number %d'%(channelNumber)
	
	configuration[channelNumber].triggerListChX = defaultTriggerListChX
	
      continue
    colonIndex = line.find(':')
    
    # If field separator not located, report syntax error
    if(colonIndex<0):
      sys.exit('Syntax error processing configuration file "%s":\n%s\n' %(configurationFile, line))
    
    parameterName = line[0:colonIndex+1]
    parameterValue = line[colonIndex+1:len(line)]
      
    colonIndex_2 = line.count(':')
    channel_flag = (parameterName, 'channelName:')
    if((colonIndex_2==1) & (channel_flag==1)):
      pv_replace_1 = parameterValue.split('_')[0]
      pv_replace_2 = pv_replace_1 + '_'
      parameterValue = parameterValue.replace(pv_replace_2, pv_replace_ + ':')
      
      parameterValue = re.sub('_[0-9]+[0-9]', '', parameterValue)
      parameterValue = re.sub('_[0-9]', '', parameterValue)
    

    parameterName = parameterName.strip()
    parameterValue = parameterValue.strip()
    
    # We assume that the first channel read from the config file is the channel H and
    # the second channel is channel X. We impose a condition to read only the first
    # channels in the configuration file.
    # Assign parameters based on name
    
    if(parameterName=='channelName:'):
      configuration[channelNumber].channelName = eval(parameterValue)
    elif(parameterName=='frameType:'):
      configuration[channelNumber].frameType = eval(parameterValue)
    elif(parameterName=='sampleFrequency:'):
      configuration[channelNumber].sampleFrequency = eval(parameterValue)
    elif(parameterName=='triggerListChH:'):
      if(channelNumber==0):
	configuration[channelNumber].triggerListChH = eval(parameterValue)
    elif(parameterName=='triggerListChX:'):
      if(channelNumber==1):
	configuration[channelNumber].triggerListChX = eval(parameterValue)
    elif(parameterName=='transferFunctionXtoH:'):
      if(channelNumber==1):
	configuration[channelNumber].transferFunctionXtoH = eval(parameterValue)
    else:
      sys.exit('Unknow configuration parameter %s %s'%(parameterName, parameterValue))
      
  
sectionIndex = np.asarray(sectionIndex)
sectionName = np.asarray(sectionName)
sectionStart = np.asarray(sectionStart)

configurationFID.close()
  
  
# Create/read channel names
  
channelHIndex = 0
channelXIndex = 1
channelYIndex = 2

# Isolate the channel number index for channels H, X and Y  
chanHName = configuration[channelHIndex].channelName
frameTypeH = configuration[channelHIndex].frameType
samplFreqH = configuration[channelHIndex].sampleFrequency

chanXName = []
samplFreqX = []
frameTypeX = []
chanXName.append(configuration[channelXIndex].channelName)
samplFreqX.append(configuration[channelXIndex].sampleFrequency)
frameTypeX.append(configuration[channelXIndex].frameType)

if(couplingModel=='bilinear'):
  chanXName.append(configuration[channelYIndex].channelName)
  samplFreqX.append(configuration[channelYIndex].sampleFrequency)
  frameTypeX.append(configuration[channelYIndex].frameType)
  # If there are two orthogonal slow channles specified in one name
  # string such as H1:ASC-QPDX_{P, Y} separate them into two strings.
  # we currently assume that there are only two orthogonal slow channels
  # This can be generalized, if needed
  startIdx = chanXName[1].find('(')
  endIdx = chanXName[1].find(')')
  midIdx = chanXName[1].find(',')
  if(startIdx>0 & midIdx >0 & endIdx>0):
    complxChannel = configuration[channelYIndex].channelName
    chanXName[1] = complxChannel[0:startIdx] + complxChannel[startIdx+1:midIdx]
    chanXName[2] = complxChannel[0:startIdx] + complxChannel[midIdx+1:endIdx]
    samplFreqX.append(samplFreqX[1])
    frameTypeX.append(frameTypeX[1])


chanXName = np.asarray(chanXName)
samplFreqX = np.asarray(samplFreqX)
frameTypeX = np.asarray(frameTypeX)
# Create output directory
outDirList = []
# Report status
print 'Creating event directory/directories...\n'
analysisStartTime = min(segStartTimeVec)
analysisEndTime  = max(segEndTimeVec)

startIdx = outDir.find('(')
endIdx = outDir.find(')')
midIdx = outDir.find(',')

if(startIdx>=0 & endIdx>=0 & midIdx>=0):
  outDirList.append(outDir[0:startIdx] + outDir[startIdx+1:midIdx])
  outDirList.append(outDir[0:startIdx] + outDir[midIdx+1:endIdx])
else:
  outDirList.append(outDir)
  
for iDir in xrange(len(outDirList)):
  os.system('mkdir -p %s'%(outDirList[iDir]))
  os.system('cp %s %s/configuration.txt'%(configurationFile, outDirList[iDir]))
  
if(debugLevel>=2):
  os.system('mkdir -p %s/debug_plots'%(outDirList[0]))

logFile = outDirList[0] + '/' + logFileName
logFid = open(logFile, 'w+')

# Write Log File Header Information
# Report Status

if(debugLevel>=0):
  logFid.write('## Bilinear coupling veto analysis [%d - %d]\n'%(analysisStartTime, analysisEndTime))
  logFid.write('## Created by %s on %s at %s\n' %(os.getenv('USER'), strftime("%Y-%m-%d"), strftime("%H:%M:%S")))
  

numTrigsH = 0
numTrigsX = 0

# Add code here if linear coupling model is to be implemented
transFnXtoH = []
if (couplingModel == 'linear'):
  transferFunctionXtoH = configuration[channelXIndex].transferFunctionXtoH
  
  if(transferFunctionXtoH != 'null'):
    transFnXtoHRead = np.loadtxt(transferFunctionXtoH)
    transFnXtoH = TransferFunctionXtoH(transFnXtoHRead[:, 0], transFnXtoHRead[:, 1])
  else:
    freq = np.linspace(0.1, 4000.0, num=4001, endpoint=True)
    transFnXtoH = TransferFunctionXtoH(freq, np.ones(len(freq), dtype=float))
    
    
    


# Read Trigger Lists for Channels H and X
triggerListChH = configuration[channelHIndex].triggerListChH

# Report status
logFid.write('LOG: Reading channel H trigger list %s...\n' %(triggerListChH))

# Open channel H trigger list file for reading
triggerListChHFID = is_non_zero_file(triggerListChH)

if(triggerListChHFID):
  # Selct only triggers passing the SNR threshold
  # Force as 2-d array
  trigDataMatrixH = np.loadtxt(triggerListChH, dtype=np.float64).reshape(-1, 5)
  triggerSignificanceH = trigDataMatrixH[:, 4]
  
  startTimeH = trigDataMatrixH[:, 0]
  endTimeH = trigDataMatrixH[:, 1]
  
  
  # Edit: No need to filter SNR threshold as the files in KWtrigs already come with
  # thresholded triggers <2018-03-30>
  trigIdxH = np.intersect1d(np.where(endTimeH <= analysisEndTime)[0],
			    np.where(startTimeH>=analysisStartTime)[0])
  # Filter out triggers whose duration is greater than 1 second
  #trigIdxH = np.intersect1d(trigIdxH, np.where((endTimeH - startTimeH)<1.0)[0])
  
  triggerListH = trigstruct.TrigStruct(trigDataMatrixH[trigIdxH, 0], trigDataMatrixH[trigIdxH, 1], trigDataMatrixH[trigIdxH, 2], trigDataMatrixH[trigIdxH, 3], triggerSignificanceH[trigIdxH])
  
  del trigDataMatrixH, triggerSignificanceH
else:
  print '### Unable to open the trigger file for channel H\n(%s)' %(triggerListChH)
  sys.exit()

# Read the trigger list for channel X
triggerListChX = configuration[channelXIndex].triggerListChX

# Report status
logFid.write('LOG: Reading channel X trigger list %s...\n' %(triggerListChX))

# Open channel X trigger list for reading
triggerListChXFID = is_non_zero_file(triggerListChX)

# Create a counter for the number of lines in channel X trigger list file
numTriggerXLines = 0

if(triggerListChXFID):
  # Force as 2-d array
  trigDataMatrixX = np.loadtxt(triggerListChX, dtype=np.float64).reshape(-1, 5)
  
  # select only triggers passing the SNR threshold
  triggerSignificanceX = trigDataMatrixX[:, 4]
  startTimeX = trigDataMatrixX[:, 0]
  endTimeX = trigDataMatrixX[:, 1]
  
  # Edit: No need to filter SNR threshold as the files in KWtrigs already come with
  # thresholded triggers <2018-03-30>
  trigIdxX = np.intersect1d(np.where(endTimeX <= analysisEndTime)[0],
			    np.where(startTimeX>=analysisStartTime)[0])
  # Filter out triggers whose duration is greater than 1 second
  #trigidxX = np.intersect1d(trigIdxX, np.where((endTimeX - startTimeX)<1.0)[0])
  
  triggerListX = trigstruct.TrigStruct(trigDataMatrixX[trigIdxX, 0],
				       trigDataMatrixX[trigIdxX, 1], trigDataMatrixX[trigIdxX, 2], trigDataMatrixX[trigIdxX, 3],
				       triggerSignificanceX[trigIdxX])
  
  del trigDataMatrixX, triggerSignificanceX
else:
  print '### Unable to open trigger file for channel X\n(%s)'%(triggerListChX)
  sys.exit(0)

nSeg = len(segStartTimeVec)

if(nSeg!=len(segEndTimeVec)):
  sys.exit('Segment start and end time vectors have diff. length. Check seg file\n')

for iSeg in xrange(nSeg):
	segStartTime = segStartTimeVec[iSeg]
  segEndTime = segEndTimeVec[iSeg]
  
  logFid.write('# Processing segment # %d .[%d, %d]\n' %(iSeg, segStartTime, segEndTime))
  
  # Select triggers in the list falling in this segment and passing significance threshold
  triggerIndexChH = np.intersect1d(np.where(triggerListH.centralTime >= segStartTime)[0], 
				   np.where(triggerListH.centralTime <= segEndTime)[0])
  
  triggerListHSeg = trigstruct.TrigStruct(triggerListH.startTime[triggerIndexChH],
					  triggerListH.endTime[triggerIndexChH],
					  triggerListH.centralTime[triggerIndexChH],
					  triggerListH.centralFrequency[triggerIndexChH],
					  triggerListH.triggerSignificance[triggerIndexChH]
					  )
  del triggerIndexChH
  
  # Find gpsStart, gspEnd and gpsCentral triggers from channel X that are not b/w
  # segStartTime and segEndTime.
  
  triggerIndexChX = np.intersect1d(np.where(triggerListX.centralTime >= segStartTime)[0], 
				   np.where(triggerListX.centralTime <= segEndTime)[0])
  
  triggerListXSeg = trigstruct.TrigStruct(triggerListX.startTime[triggerIndexChX],
					  triggerListX.endTime[triggerIndexChX],
					  triggerListX.centralTime[triggerIndexChX],
					  triggerListX.centralFrequency[triggerIndexChX],
					  triggerListX.triggerSignificance[triggerIndexChX]
					  )
  del triggerIndexChX
  
  numTrigsHseg = len(triggerListHSeg.centralTime)
  numTrigsXseg = len(triggerListXSeg.centralTime)
  
  numTrigsH = numTrigsH + numTrigsHseg
  numTrigsX = numTrigsX + numTrigsXseg
  
  # Perform Timeshift and veto analysis
  
  if((numTrigsHseg > 0) & (numTrigsXseg >0)):
    # Construct a time shift vector, make sure that the zero lag is included
    # perform the veto analysis for each time shift (incl. zero lag)
    timeShiftVec = np.linspace(timeShiftMin, timeShiftMax, numTimeShifts)
    timeShiftVec = np.unique(np.append([0], np.round(timeShiftVec)))
    
		COINC_TIME_WINDOW = 1.0

		segLength = 3600

		uniqueArgument = 'nonunique'
    
		[coincTrigH, coincTrigX] = bcv.mcoinc(maxNumCoinc, triggerListHSeg.centralTime, triggerListXSeg.centralTime, COINC_TIME_WINDOW, segLength, uniqueArgument)
		if(len(coincTrigH)==len(coincTrigX) & len(coincTrigH)>0):
			#Read data for timeShift from frameCache
			logFid.write('LOG: Number of coincidences : %d...\n' %(len(coincTrigH)))
			
			for iTimeShift in xrange(len(timeShiftVec)):
				timeShift = timeShiftVec[iTimeShift]
				vetoanalysis.vetoanalysis(frameCacheFileH, frameCacheFileX, [chanHName], chanXName, [frameTypeH], frameTypeX, samplFreqH, samplFreqX, coincTrigH, coincTrigX,
				highPassCutoff, triggerListHSeg, triggerListXSeg,
				couplingModel, transFnXtoH, segStartTime, segEndTime,
				timeShift, outDirList, logFid, debugLevel)
		else:
			logFid.write('LOG: No coincident triggers found in this segment')
  else:
    logFid.write('LOG: No triggers in this segment numTrigsHseg = %d numTrigsXseg = %d\n'%(numTrigsHseg, numTrigsXseg))
  
  del triggerListHSeg, triggerListXSeg

for iDir in xrange(len(outDirList)):
  textSummaryFID = open(outDirList[iDir] + '/summary.txt', 'w')
  
  # print summary of results
  textSummaryFID.write('# Summary file of veto analysis')
  textSummaryFID.write('Created by : %s\n' %(os.getenv('USER')))
  textSummaryFID.write('Created on: %s\n'%(strftime("%Y-%m-%d")))
  textSummaryFID.write('analysisStartTime UTC : %s\n' %(gtime.from_gps(analysisStartTime)))
  textSummaryFID.write('segEndTime UTC : %s\n' %(gtime.from_gps(segEndTime)))
  textSummaryFID.write('configurationFile : %s\n' %(configurationFile))
  textSummaryFID.write('frameCacheFileH : %s\n'%(frameCacheFileH))
  textSummaryFID.write('frameCacheFileX : %s\n'%(frameCacheFileX))
  textSummaryFID.write('couplingModel :%s\n' %(couplingModel))
  textSummaryFID.write('highPassCutoff : %3.2f\n' %(highPassCutoff))
  textSummaryFID.write('outDir : %s\n' %(outDirList[iDir]))
  textSummaryFID.write('logFile : %s\n' %(logFile))
  textSummaryFID.write('debugLevel : %d\n' %(debugLevel))
  textSummaryFID.write('numTrigsH :%d\n' %(numTrigsH))
  textSummaryFID.write('numTrigsX: %d\n' %(numTrigsX))
  
  textSummaryFID.close()

# Exit

# Report completion
if(debugLevel>=0):
  logFid.write('Finished on %s at %s\n'%(strftime("%Y-%m-%d"), strftime("%H:%M:%S")))

logFid.close()

for iDir in range(1, len(outDirList)):
  os.system('cp %s %s/%s'%(logFile, outDirList[iDir], logFileName))


      
      
   
    
  

  


  
    
  
	
    
    
    
    
    
    
    
    
    
    
    





  
