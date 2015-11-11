# Bilinear coupling veto (bcv) module
# this module contains functions to read, write and process data obtained from the 
# detector channel and auxilary channels.
# The functions are called by the main script (vetoanalysis.py) to process this data
# based on inputs received by the Event Trigger Generator (ETG) Algorithms

# The firls() function has been retrieved from python ticket 648
# http://projects.scipy.org/scipy/attachment/ticket/648/designtools.py
# 
# The multirate code has been used from the following site
# http://mubeta06.github.io/python/sp/_modules/sp/multirate.html

from pylal import Fr
import sys
import scipy.signal as sig
import numpy as np
import scipy.linalg as linalg
import scipy.interpolate as sinterp

def readData(frameCache, channelNames,frameTypes, startTime, stopTime, timeShifts, debugLevel):
  numberOfChannels = len(channelNames)
  if(len(frameTypes)!=numberOfChannels):
    print channelNames, frameTypes
    sys.exit('Number of frame types is inconsistent with number of channels')
  
  # This function reads data from the file frameCache which contains data in channelName
  
  # Permit redundant data
  allowRedundantFlag = True
  # Initializing timeShifts if they are null
  if(len(timeShifts)==0):
    timeShifts = np.zeros(len(numberOfChannels))

  if(len(timeShifts)!=numberOfChannels):
    sys.exit('Number of time shifts inconsistent with number of channels')
  # Making sure stop time greater than start time
  if(stopTime<startTime):
    stopTime = startTime + stopTime
  data = []
  samplFreq = []
  
  #Reading data for all channels one by one and storing them in vectors
  for channelNumber in xrange(numberOfChannels):
    # The function frgetvect is from the python implementation of LALSuite and returns
    # the following
    # frame[0] is an array containing the strains
    # frame[1] is the end time of the data
    # frame[3] is the inverse of the Sampling Frequency
    #print "Reading data of channel %s\n" %(channelNames[channelNumber])
    mData, mSamplFreq = readframedata(frameCache, channelNames[channelNumber],
				       frameTypes[channelNumber], startTime - timeShifts[channelNumber], stopTime - timeShifts[channelNumber],
				       allowRedundantFlag, debugLevel)
    data.append(mData)
    samplFreq.append(mSamplFreq)
  data = np.asarray(data)
  samplFreq = np.asarray(samplFreq)
  return [data, samplFreq]



#def resample2(data, samplFreqD, samplFreq):
  #dataDuration = len(data)/samplFreqD
  #from fractions import Fraction
  #import sp_multirate as sp
  ## (http://mubeta06.github.io/python/sp/_modules/sp/multirate.html)
  
  #frac = Fraction(samplFreq/samplFreqD).limit_denominator(1000)
  #upSampleFactor = frac.numerator
  #downSampleFactor = frac.denominator
  #filterOrder = 2*256*max(upSampleFactor, downSampleFactor)
  #filterCutOff = 0.99/max(upSampleFactor, downSampleFactor)
  #filterFrequencies = np.asarray([0, filterCutOff, filterCutOff, 1])
  ##print 'filterFrequencies: ', filterFrequencies
  #filterMagnitudes = np.asarray([1, 1, 0, 0])
  #filterCoefficients = upSampleFactor*firls(filterOrder, filterFrequencies, filterMagnitudes)*sig.hanning(filterOrder+1)
  ##print 'filterCoefficients: ', filterCoefficients
  
  #dataInter = sp.resample(data, upSampleFactor, downSampleFactor, filterCoefficients)
  ##print 'dataInter: ', dataInter
  #return dataInter

def resample2(data, samplFreqD, samplFreq):
  totalTime = len(data)/samplFreqD
  new_num = totalTime*samplFreq
  dataInter = sig.resample(data, new_num)
  return dataInter

def linearCouplingCoeff(dataH, dataX, timeH, timeX, transFnXtoH, segStartTime,
			segEndTime, timeShift, samplFreq, logFid, debugLevel):
  # LINEARCOUPLINGCOEFF - calculate the cross correlation coeff b/w the gravitational
  # ave channel H and the "projected" instrumental channel X. The noise in the
  # instrumental channel X is projected to the domain of the H using a linear coupling
  # function Txh

  MIN_FREQ = 10.0
  MAX_FREQ = 4000.0  
  IFO_LENGTH = 4000
  
  if((len(dataH)==0) | (len(dataX)==0)):
    logFid.write('Error: One or more data vectors are empty..\n')
    logFid.write('Error: len(dataH) = %d len(dataX) = %d..\n' %(len(dataH), len(dataX[0])))
  
  elif(len(dataH)!=len(dataX[0])):
    logFid.write('Error: Different lengths. len(dataH) = %d len(dataX) = %d..\n'%(len(dataH), len(dataX[0])))
  else:
    dataH = dataH - np.mean(dataH)
    dataX = dataX[0] - np.mean(dataX[0])
    
    segIdxH = np.intersect1d(np.where(timeH>=segStartTime)[0], np.where(timeH<segEndTime)[0])
    dataH = dataH[segIdxH]
    
    segIdxX = np.intersect1d(np.where(timeX + timeShift >= segStartTime)[0], np.where(timeX + timeShift < segEndTime)[0])
    dataX = dataX[segIdxX]
    
    nfft = len(dataH)
    
    freqVecH = np.fft.rfftfreq(nfft, 1.0/samplFreq)
    fftChanH = np.fft.rfft(dataH)
    
    nfft = len(dataX)
    freqVecX = np.fft.rfftfreq(nfft, 1.0/samplFreq)
    fftChanX = np.fft.rfft(dataX)
    
    freqBandIdx = np.intersect1d(np.where(transFnXtoH.frequency>=MIN_FREQ)[0], np.where(transFnXtoH.frequency < MAX_FREQ)[0])
    if(len(freqBandIdx)!=0):
      transFnXtoH.frequency = transFnXtoH.frequency[freqBandIdx]
      transFnXtoH.Txh = transFnXtoH.Txh[freqBandIdx]
    
    TxhFreqMin = np.min(transFnXtoH.frequency)
    TxhFreqMax = np.max(transFnXtoH.frequency)
    
    
    #if(len(freqBandIdx)!=0):
      #transFnXtoH.frequency = transFnXtoH.frequency[freqBandIdx]
      #transFnXtoH.Txh = transFnXtoH.Txh[freqBandIdx]
    if(len(freqVecH)!=len(freqVecX)):
      logFid.write('ERROR: Unequal size for freq. vectors.\n')
      logFid.write('ERROR: len(freqVecH) = %d len(freqVecX) = %d.\n' %(length(freqVecH), length(freqVecX)))
    else:
      freqBandIdx = np.intersect1d(np.where(freqVecH>=TxhFreqMin)[0], np.where(freqVecH < TxhFreqMax)[0])
      
      if(len(freqBandIdx)!=0):
	fftChanH = fftChanH[freqBandIdx]
	fftChanX = fftChanX[freqBandIdx]
	freqVecH = freqVecH[freqBandIdx]
	freqVecX = freqVecX[freqBandIdx]
      
      
      freqResolTxh = np.float(samplFreq)/len(dataX[0])
      print 'freqResolTxh ', freqResolTxh
      
      [tFreqIntp, tfMagIntp,tfPhaseIntp] = interpolatetransfn(transFnXtoH.frequency,
							     np.abs(transFnXtoH.Txh),
							     np.unwrap(np.angle(transFnXtoH.Txh)), freqResolTxh)
      TxhInterp = tfMagIntp* np.exp(1j*tPhaseIntp)
      
      if(len(fftChanX)==size(TxhInterp)):
	xPrime = fftChanX*TxhInterp/IFO_LENGTH
      else:
	logFid.write('ERROR: size(fftChanX) = %d size(TxhInterp) = %d\n' %(len(fftChanX), len(TxhInterp)))
	logFid.write('ERROR: fftChanX and TxhInterp have different sizes.\n')
	sys.exit('Inconsistent dimensions of data and transfer function')
      
      [rXH, rMaxXH] = calcrossCorr(xPrime, fftChan)
      return [rH, rMaxXH]
    
    
  
  
def bilinearCouplingCoeff(dataH, dataP, timeH, timeP,
			  segStartTime,segEndTime, timeShift, samplFreq, logFid, debugLevel):
  # This function computes the correlation coefficient between Channel H and al
  # the pseudo channels P.
  #
  # Usage: [rPH, rPHAbs] = bilinearCouplingCoeff(dataH, dataP, timeH, timeP,
  #                               segStartTime, segEndTime, timeShift, samplFreq
  #                               logFid, debugLevel)
  # Set the frequency range of the veto analysis
  
  
  MIN_FREQ = 10.0
  MAX_FREQ = 4000.0
  
  # Meta Data
  segIdxH = np.intersect1d(np.where(timeH>=segStartTime)[0], np.where(timeH<segEndTime)[0])
  dataH = dataH[segIdxH]
  dataH = dataH - np.mean(dataH)

  
  # Set parameters for calculating spectrogram
  nfft = len(dataH)
  #print "nfft = %d" %(nfft)
  wind = np.ones(nfft)
  
  # Calculate spectrogram
  #print 'dataH: ', dataH
  #print 'len(dataH): ', len(dataH)
  freqVecH = np.fft.rfftfreq(nfft, 1.0/samplFreq)
  fftChanH = np.fft.rfft(dataH)
  #print 'fftChanH: ', fftChanH
  #print 'len(fftChanH): ', len(fftChanH)


  [numChanP, lengthP] = np.shape(dataP)
  segIdxP = np.intersect1d(np.where(timeP + timeShift >= segStartTime)[0], np.where(timeP + timeShift < segEndTime)[0])
  rPH = np.asarray([])
  rPHAbs = np.asarray([])
  
  #print 'numChanP: ', numChanP
  for iChan in xrange(numChanP):
    dataVecP = dataP[iChan]
    dataVecP = dataVecP[segIdxP]
    dataVecP = dataVecP - np.mean(dataVecP)
    
    freqVecP = np.fft.rfftfreq(nfft, 1.0/samplFreq)
    fftChanP = np.fft.rfft(dataVecP)
    #print 'fftChanP.shape: ', fftChanP.shape, 'fftChanH.shape', fftChanH.shape
    
    # Make sure the fft vectors are of the same size
    if(len(freqVecP)!=len(freqVecH)):
      logFid.write('ERROR: Unequal size of freq. vectors. \n')
      logFid.write('ERROR: len(freqVecH) = %d len(freqVecP) = %d\n', len(freqVecH), len(freqVecP))
      return
    else:
      #Select the frequency band
      freqBandIdx = np.intersect1d(np.where(freqVecH>MIN_FREQ)[0], np.where(freqVecP < MAX_FREQ)[0])
      
      # Calculate the cross-correlation statistic for the segment of data in channels h
      # and "projected" P
      #print "freqVecP = ", freqVecP
      #print "fftChanP.shape", fftChanP.shape
      [a, b] = calCrossCorr(fftChanP[freqBandIdx], fftChanH[freqBandIdx])
      #print 'a: ', a, 'b: ', b
      rPH = np.append(rPH, a)
      rPHAbs = np.append(rPHAbs, b)
  
  return [rPH, rPHAbs]

def interpolatetransfn(tfFreq, tfMag, tfPhase, reqFreqRes):
  
  # Compute the required frequency resolution.
  tfFRes = tfFreq[1] - tfFreq[0]
  print 'tfFRes ', tfFRes
  
  # Create a new frequency vector for interpolating the transfer function.
  newFreqVec = np.arange(np.min(tfFreq), np.max(tfFreq)+ reqFreqRes, reqFreqRes)
  #Use the resample command for interpolating the transfer function.
  #tfFreqIntp  = newFreqVec
  #tfMagIntp   = resample(tfMag,tfFRes,reqFreqRes)
  #tfPhaseIntp = resample(tfPhase,tfFRes,reqFreqRes)
  tck = sinterp.splprep(tfFreq, tfFreq, s=0)
  tfFreqIntp = sinterp.splev(newFreqVec, tck, der=0)
  
  tck = sinterp.splprep(tfFreq, tfMag, s=0)
  tfMagIntp = sinterp.splev(newFreqVec, tck, der=0)
  
  tck = sinterp.splprep(tfFreq, tfPhase, s=0)
  tfPhaseIntp = sinterp.splev(newFreqVec, tck, der=0)
  
  return [tfFreqIntp, tfMagIntp,tfPhaseIntp]
  

def calCrossCorr(u,v):
  # Auxlilary function which calculates the cros correlation between two signals
  # The inputs are the spectrograms of the two signals
  if(u.shape!=v.shape):
    sys.exit('Size of u and v mist be the same\n')
      
  uDotu = np.sum(u*np.conj(u))
  vDotv = np.sum(v*np.conj(v))
  uConjv = u*np.conj(v)/(np.sqrt(uDotu)*np.sqrt(vDotv))
  #print 'uDotu :', uDotu, 'vDotv: ', vDotv, 'u*u/v*v :', uDotu/vDotv
  r = np.real(np.sum(uConjv))
  
  xCorr = np.fft.irfft(uConjv)/len(u)
  rMax = np.max(np.real(xCorr))
  rMin = np.min(np.real(xCorr))
  
  rAbs = rMax
  if(np.abs(rMin)>np.abs(rMin)):
    rAbs = np.abs(rMin)
  return [r, rAbs]

def highpass(rawData, samplFreq, highpassCutOff):
  #This function high passes the data with sampling frequency samplFreq with a highpass
  # cut off of highpassCutOff
  #
  # Usage: [highPassedData] = highpass(rawData, samplFreq, highpassCutOff)
  # 
  # rawData : raw time series Data
  # samplFreq: sampling Frequency of Data
  # highpassCutOff: The high pass frequency cut off.
  #
  numberOfChannels = len(rawData)
  dataLength = len(rawData[0])
  duration = dataLength/samplFreq
  
  
  halfDataLength = dataLength/2 + 1
  
  for channelNumber in xrange(numberOfChannels):
    if(len(rawData[channelNumber]) != dataLength):
      sys.exit('Data length not consistent\n')
  
  nyquistFrequency = samplFreq/2.0
  
  lpefOrder  = 0
  
  if(highpassCutOff>0):
    hpfOrder = 12
    hpfZeros, hpfPoles, hpfGain = sig.butter(hpfOrder, highpassCutOff/nyquistFrequency, btype = 'highpass', output = 'zpk' )
    hpfSOS = sig.zpk2sos(hpfZeros, hpfPoles, hpfGain)
    
    #magnitude response of high pass filter
    minimumFrequencyStep = 1.0/duration
    frequencies = np.arange(0, nyquistFrequency, minimumFrequencyStep)
    hpfArgument = np.power((frequencies / highpassCutOff), 2*hpfOrder)
    hpfResponse = hpfArgument/(1 + hpfArgument)
    
    highPassCutOffIndex = np.ceil(highpassCutOff/minimumFrequencyStep)
  
  highPassedData = []
  for channelNumber in xrange(numberOfChannels):
    if(highpassCutOff>0):
      x = sig.sosfilt(hpfSOS, rawData[channelNumber])
      x = np.flipud(x)
      x = sig.sosfilt(hpfSOS, x)
      x = np.flipud(x)
    else:
      x = rawData[channelNumber]
    
    x[0:lpefOrder] = np.zeros(lpefOrder)
    x[dataLength - lpefOrder:dataLength-1] = np.zeros(lpefOrder)
    
    highPassedData.append(x)
  
  highPassedData = np.asarray(highPassedData)
  
  return highPassedData
  
class FrameCacheStruct:
  def __init__(self, sites, frameTypes, startTimes, stopTimes, durations, directories ):
    self.sites = sites
    self.frameTypes = frameTypes
    self.startTimes = startTimes
    self.stopTimes  = stopTimes
    self.durations  = durations
    self.directories = directories

def loadframecache(cachePath):
  #loadframecache loads frame file cache information
  
  #This function reads the lal/frame file cache information stored in the 
  #specified file. The resulting cache structure is used to locate frame data
  #during subsequent calls to readframedata
  
  #usage: cache = loadframecache(cachePath)
  
  #cachePath: path to lal/frame cachef file
  #cache: lal/frame cache structure
  
  #formats:
    
  #FRAMECACHE: The frame cache file should consist of whitespace delimited ASCII
  #text and contains one line for each contiguous data segment with a common site, frame
  #type, duration and directory. Each line should consist of the following six columns
  
    #* site designator (e.g 'H' or 'L')
    #* frame file type (e.g 'RDS_R_L3')
    #* GPS start time of segment
    #* GPS stop time of segment
    #* frame file duration in seconds
    #* full path name of directory
  #The resulting cache frame structure consists of the following six fields.
  
    #* .sites              numpy array of segment site designators
    #* .frameTypes         numpy array of segment frame file frameTypes
    #* .startTimes         numpy array of segment GPS start times
    #* .stopTimes          numpy array of segment GPS stop times
    #* .durations          numpy array of segment frame file durations
    #* .directories        numpy array of segment directory names
    
  sites     = np.genfromtxt(cachePath,delimiter=' ',usecols=(0), dtype=None,unpack=True)
  frameTypes = np.genfromtxt(cachePath,delimiter=' ',usecols=(1), dtype=None,unpack=True)
  startTimes = np.genfromtxt(cachePath,delimiter=' ',usecols=(2), dtype=None,unpack=True)
  stopTimes = np.genfromtxt(cachePath,delimiter=' ',usecols=(3), dtype=None,unpack=True)
  durations = np.genfromtxt(cachePath,delimiter=' ',usecols=(4), dtype=None,unpack=True)
  directories =np.genfromtxt(cachePath,delimiter=' ',usecols=(5), dtype=None,unpack=True)
  
  cache = FrameCacheStruct(sites, frameTypes, startTimes, stopTimes,durations, directories)
  
  return cache
  
  
def  readframedata(frameCache, channelName, frameType, startTime, stopTime,
		   allowRedundantFlag, debugLevel):
  #READFRAMEDATA Read a single channel of data from frame files
  #
  #READFRAMEDATA finds and retrieves the requested time series data from a
  #set of frame files.  The data is specified by the frame file type,
  #channel name, start time, and duration or stop time.  The necessary frame
  #files are located using a file name caching scheme.
  #%
  #usage: [data, sampleFrequency, time] = ...
  #readframedata(frameCache, channelName, frameType, ...
  #startTime, stopTime, allowRedundantFlag, ...
  #debugLevel);
  #%
  #frameCache           file name cache
  #channelName          channel name
  #frameType            frame file type
  #startTime            GPS start time
  #stopTime             GPS stop time (or duration)
  #allowRedundantFlag   permit redundant frame data
  #debugLevel           verboseness of debug output
  #%
  #data                 data vector
  #sampleFrequency      sample frequency [Hz]
  #time                 time vector
  #%
  #READFRAMEDATA expects frame cache information in the format produced by
  #LOADFRAMECACHE which contains the site, frame type, duration, time, and
  #location of available frame data files.
  #%
  #The requested site designator is determined from the first character of
  #the requested channel name.  The frame file type may contain wildcards.
  #Unless redundant frame data is permitted, it is an error if the same
  #frame file appears more than once in the frame cache file.  If redundant
  #frame data is permitted, the first matching frame file from the cache is
  #used.  It is always an error if two frame files overlap but do not cover
  #the same time range or differ in type.  By default, redundant frame data
  #is permitted.
  #%
  #READFRAMEDATA retrieves data from the requested start time up to, but not
  #including, the requested stop time, such that stop minus start seconds
  #are retrieved.  Alternatively, the desired duration in seconds may be
  #specified instead of the GPS stop time parameter.
  #%
  #The resulting time series is returned as two row vectors containing the
  #data sequence and the corresponding GPS timestamps, as well as the scalar
  #sample frequency.  To protect against roundoff error, an integer sample
  #frequency is assumed.
  #%
  #If it is unable to load the requested data, READFRAMEDATA returns empty
  #result vectors and zero sample frequency as well as a warning if
  #debugLevel is set to 1 or higher.  By default, a debugLevel of unity is
  #assumed.
  #%
  #READFRAMEDATA is built on top of the FRGETVECT function from the FrameL
  #library, which is available from the following URL.
  #%
  #
  #%
  #
  #
  #Shourov K. Chatterji <shourov@ligo.mit.edu>
  #Jameson Rollins <jrollins@phys.columbia.edu>
  #
  #$Id: readframedata.m 2326 2009-09-21 08:37:42Z jrollins $
  #
  # Rewritten in Python by Sudarshan Ghonge <sudu.ghonge@gmail.com>
  # 2015-09-04

  #
  # if specified stop time precedes start time,
  if(stopTime<startTime):
    # treat stop time as a duration
    stopTime = startTime + stopTime
  
  # determine site designator from channel name
  site = channelName[0]
  
  #if specified frame cache is invalid
  if(frameCache==None):
    if(debugLevel>=1):
      # Issue warning
      print 'Warning: Invalid frame cache'
    
    data = []
    time = []
    sampleFrequency = 0
    return [data, sampleFrequency]
    # return empty results
    
  
  # Identifying matching segments from frame cache
  
  # find overlap of cache segments with requested data
  segmentStartTimes = np.maximum(startTime, frameCache.startTimes)
  segmentStopTimes = np.minimum(stopTime, frameCache.stopTimes)
  
  # identify cache segments which overlap requested times
  segments = np.where(segmentStopTimes > segmentStartTimes)[0]
  
  # if no segments overlap with requested times
  if(len(segments)==0):
    
    if debugLevel>=1:
      print 'Warning: No data available for [%d, %d]' %(startTime, stopTime)
    data = []
    time = []
    sampleFrequency = 0
    # return empty results
    return [data, sampleFrequency]
  # otherwise, find overlapping segments
  else:
    # identify cache segments with requested site and frame type
    siteMatches = []
    sitesScanned=0
    for iSite in frameCache.sites[segments]:
      if site in iSite:
	siteMatches.append(sitesScanned)
      sitesScanned+=1
      
    siteMatches = np.asarray(siteMatches)
    
    frameTypeMatches = []
    frameTypesScanned=0
    print 'frameCache.frameTypes[segments]: ', frameCache.frameTypes[segments]
    print 'frameType: ', frameType
    for iType in frameCache.frameTypes[segments]:
      if frameType in iType:
	frameTypeMatches.append(frameTypesScanned)
      frameTypesScanned+=1
    
    frameTypeMatches = np.asarray(frameTypeMatches)
    
    segIdx = np.intersect1d(siteMatches, frameTypeMatches)
    if(len(segIdx)==0):
      segments = []
    else:
      segments = segments[segIdx]
  
  # Identify available frame files
  
  # initialize list of available frame files
  frameFilePaths = []
  frameFileTypes = []
  frameFileStartTimes = []
  frameFileStopTimes = []
  
  # lopp over the matching segments
  for segment in segments:
    
    # frame type of the frame files in segment
    frameFileType = frameCache.frameTypes[segment]
    
    firstFrameFileStartTime = frameCache.startTimes[segment] + frameCache.durations[segment]*np.floor((segmentStartTimes[segment]-frameCache.startTimes[segment])/frameCache.durations[segment])
    
    
    lastFrameFileStartTime = frameCache.startTimes[segment] + frameCache.durations[segment]*np.ceil((segmentStopTimes[segment] - frameCache.startTimes[segment])/frameCache.durations[segment] - 1)
    
    for frameFileStartTime in np.arange(firstFrameFileStartTime,
					lastFrameFileStartTime + frameCache.durations[segment],
					frameCache.durations[segment]):
      
      frameFileStopTime = frameFileStartTime +  frameCache.durations[segment]
      
      frameFilePath = frameCache.directories[segment] + '/' + frameCache.sites[segment] + '-' + frameCache.frameTypes[segment] + '-' + '%09d' %(frameFileStartTime) + '-' + '%d' %(frameCache.durations[segment]) + '.gwf'
      
      import os.path
      if(not os.path.isfile(frameFilePath)):
	frameFilePath = frameCache.directories[segment] + '/' + frameCache.sites[segment] +'-' + frameCache.frameTypes[segment] + '-' + '%010d' %(frameFileStartTime) + '-' + '%d'%(frameCache.durations[segment]) + '.gwf'
      
      if(os.path.isfile(frameFilePath)):
	frameFilePaths.append(frameFilePath)
	frameFileTypes.append(frameFileType)
	frameFileStartTimes.append(frameFileStartTime)
	frameFileStopTimes.append(frameFileStopTime)
  
  frameFilePaths = np.asarray(frameFilePaths)
  frameFileTypes = np.asarray(frameFileTypes)
  frameFileStartTimes = np.asarray(frameFileStartTimes)
  frameFileStopTimes = np.asarray(frameFileStopTimes)
  
  numberOfFrameFiles = len(frameFilePaths)
  
  
  keepFrameFileNumbers = []
  
  for frameFileNumber in range(numberOfFrameFiles):
    keepFrameFileFlag = True
    
    for previousFrameFileNumber in range(frameFileNumber):
      
      overlapStartTime = np.maximum(frameFileStartTimes[frameFileNumber],
				    frameFileStartTimes[previousFrameFileNumber])
      overlapStoptime = np.minimum(frameFileStopTimes[frameFileNumber],
				   frameFileStopTimes[previousFrameFileNumber])
      
      if (overlapStartTime < overlapStoptime):
	if(allowRedundantFlag):
	  if( (frameFileStartTimes[frameFileNumber]==frameFileStartTimes[previousFrameFileNumber])      & (frameFileStopTimes[frameFileNumber] == frameFileStopTimes[previousFrameFileNumber])
          & (frameFileTypes[frameFileNumber]==frameFileTypes[previousFrameFileNumber])):
	    keepFrameFileFlag = False
	    continue
	  else:
	    if(debugLevel>=1):
	      print 'Warning: Overlapping but dissimilar frame files %s and %s.' %(frameFilePaths[frameFileNumber], frameFilePaths[previousFrameFileNumber])
	    data = []
	    time = []
	    sampleFrequency=0
	    return [data, sampleFrequency]
	else:
	  if(debugLevel>=1):
	    print 'Warning: Redundant frame files %s and %s.' %(frameFilePaths[frameFileNumber], frameFilePaths[previousFrameFileNumber])
	  data = []
	  time = []
	  sampleFrequency = 0
	  return [data, sampleFrequency]
	
    if(keepFrameFileFlag):
      keepFrameFileNumbers.append(frameFileNumber)
  
  keepFrameFileNumbers = np.asarray(keepFrameFileNumbers)
  frameFilePaths = frameFilePaths[keepFrameFileNumbers]
  frameFileTypes = frameFileTypes[keepFrameFileNumbers]
  frameFileStartTimes = frameFileStartTimes[keepFrameFileNumbers]
  frameFileStopTimes = frameFileStopTimes[keepFrameFileNumbers]
  
  sortedIndices = np.argsort(frameFileStartTimes)
  frameFilePaths = frameFilePaths[sortedIndices]
  frameFiletypes = frameFileTypes[sortedIndices]
  frameFileStartTimes = frameFileStartTimes[sortedIndices]
  frameFileStopTimes = frameFileStopTimes[sortedIndices]
  
  continuityStartTimes = np.append(np.maximum(startTime, frameFileStartTimes), stopTime)
  continuityStopTimes = np.append(startTime, np.minimum(stopTime, frameFileStopTimes))
  discontinuities =  np.where(continuityStartTimes!=continuityStopTimes)[0]
  
  if(len(discontinuities)>0):
    if(debugLevel >=1):
      print 'Warning: Missing %s '%(channelName), frameType[1:len(frameType)-1], ' data at ' , '%d'%(np.round(continuityStopTimes[discontinuities[0]])), '.'
    data = []
    time = []
    sampleFrequency = 0
    return [data, sampleFrequency]
  
  data = np.array([])
  time = np.array([])
  sampleFrequency = None
  
  numberOfFrameFiles = len(frameFilePaths)
  for frameFileNumber in range(numberOfFrameFiles):
    frameFilePath = frameFilePaths[frameFileNumber]
    
    if(debugLevel>=2):
      print 'Reading %s...\n' %(frameFilePath)
    
    frameFileStartTime = frameFileStartTimes[frameFileNumber]
    
    frameFileStopTime = frameFileStopTimes[frameFileNumber]
    
    readData = []
    
    readStartTime = np.maximum(startTime, frameFileStartTime)
    
    readDuration = np.minimum(stopTime, frameFileStopTime) - readStartTime
    
    realChannelName = channelName
    
    try:
      outputStruct = Fr.frgetvect(frameFilePath, realChannelName, readStartTime, readDuration, False)
      readData = outputStruct[0]
      readTime = outputStruct[2]
      readSampleFrequency = 1.0/outputStruct[3][0]
      readTimeStep = outputStruct[3][0]
      readGPS = outputStruct[1]
      
    except Exception as inst:
      if(debugLevel>=2):
	sys.exit( inst.message)
    if((len(readData)==0) | np.any(np.isnan(readData))):
      if(debugLevel>=1):
	print 'Warning: Error reading %s from %s.' %(channelName, frameFilePath)
	
      data = []
      time = []
      sampleFrequency = 0
      return [data, sampleFrequency]
    if(sampleFrequency==None):
      sampleFrequency = readSampleFrequency
    elif(sampleFrequency!=readSampleFrequency):
      if(debugLevel>=1):
	print 'Warning: Inconsistent sample frequency for %s in frameFilePath.' %(channelname, frameFilePath)
      data = []
      time = []
      sampleFrequency = 0
      return [data, sampleFrequency]
    
    data = np.append(data, readData)
  
  return [data, sampleFrequency]

#def firls(N, frequencies, pass):
  
  #weight = np.ones(len(pass)/2.0)
  #str = []
  #if(len(frequencies)!=len(pass)):
    #sys.exit('F and A must have equal lengths.')
  
  #N+=np.mod(N,2)
  
  #M = N/2
  #w_len = len(weight)
  #p_len = len(pass)
  #w = np.kron(weight.reshape(w_len, 1), np.array([[-1], [1]]))
  #omega = np.array([frequencies*np.pi])
  #i1 = np.arange(1:p_len+1, 2)
  #i2 = np.arange(2:p_len+1, 2)
  
  ## Generate the matrix Q
  #cos_ints = np.append(omega, np.sin(np.kron(np.arange(1, N+1).reshape(N, 1), omega)), axis=0)
  #q = np.append([1], [1.0/np.arange(1, 5)])*(np.dot(cos_ints,w).reshape(1, N+1)[0])
  #Q = linalg.toeplitz(q[0:M+1] + linalg.hankel(q[0:M+1], q[M:len(q)]))
  
  #omega = omega[0]
  #cos_ints2 = omega[i1]**2.0 - omega[i2]**2.0
  

# Python module which contains the following functions
# mcoinc, mhighpass mSpecgram

def mcoinc(maxCoinc, A, B, P, seglen, *args):
  # This function looks for coincident triggers b/w X and H channels
  # Assumes first row of A and B are time vectors
  # The matrices are then split in to submatrices of length seglen for processing.
  # This results in large speed
  # usage : [aidx, bidx] = mcoinc(maxcoinc, A, B, P, seglen)
  
  # A - Vector containing the central time co-ordinates of the triggers in H channel
  # B - Vector containing the central time co-ordinates of the triggers in X channel
  # P - Scalar describing the coincidence time window vector
  #
  # Example
  #
  # [hldx, xldx] = mcoinc(maxcoinc, hTimeVec, xTimeVec, timeWind, 3600, 'nonunique')
  #
  # In the above example, xTimeVec(xldx) will give the triggers in X that are coincident with the triggers in H.
  # timeWind is the time window usually, 0.5 seconds
  # maxCoinc is the max no of coincidences that are allowed (usually for logistical reasons)
  # such as memory limitations
  # If there are no such limitions, set it to 
  # maxCoinc = max(length(xTimeVec, length(hTimeVec)))
  #
  # Sudarshan Ghonge <sudu.ghonge@gmail.com>
  # Based on code written by
  # M. Hewitson
  # Aaron B. Pearlman <aaronpl@umbc.edu>
  # P. Ajith
  
  
  
  # Find the minimum of time contained in the first row of the matrices A and B
  t0 = np.int(np.floor(min(np.min(A), np.min(B) ) ))
  # Find the max of time contained in the first row of the matrices A and B
  tmax = np.int(np.round(max(np.max(A), np.max(B))))
  
  # Initialize the output row vectors that will indicate coincident trigger
  # columns b/w matrices A and B
  c1 = np.empty([0, 1], dtype=int)
  c2 = np.empty([0, 1], dtype=int)
  
  # Uniqueness argument defaulted to False
  nu = False
  
  if(len(args)>0):
    if(args[0]=='nonunique'):
      nu = True
    elif(args[0]=='unique'):
      nu = False
    else:
      sys.exit('Invalid option for variable argument in mocoinc(maxCoinc, A, B, P, seglen, *args)')
  
  for ts in xrange(t0, tmax, seglen):
    
    # get the index of this segment of the A matrix
    idxa = np.intersect1d(np.where(A >= ts)[0], np.where(A < ts + seglen)[0])
    ATemp = A[idxa]
    # get the index for this segment of the B matrix
    idxb = np.intersect1d(np.where(B >= ts)[0], np.where(B< ts + seglen)[0])
    BTemp = B[idxb]
    
    # Use the mfindcoinc function to identify coincident trigger columns
    # between matrices ATemp and BTemp
    [C1, C2] = mfindcoinc(maxCoinc, ATemp, BTemp, P)
    
    if(len(C1)>0):
      if nu:
	c1 = np.append(c1, C1 + idxa[0])
      else:
	c1 = np.append(c1, np.unique(C1) + idxa[0])
    if(len(C1)>0):
      if nu:
	c2 = np.append(c2, C2 + idxb[0])
      else:
	c2 = np.append(c2, np.unique(C2) + idxb[0])
  
  return [c1, c2]
	
	
def mfindcoinc(maxCoinc, A, B, P):
  # Auxilary function used by mcoinc() to find co-incidences between time vectors
  # A and B which are usually segments of the larger H and X channel trigger time
  # co-ordniates
  #
  # Inputs: 
  # maxCoinc - the maximum number of co-incidences
  # A  - segment of the vector containing the time coordinates (seconds) triggers from a channel
  # B  - same as above
  # P - size of the coincidence window
  
  n=0
  cont=True
  C1 = []
  C2 = []
  # Iterate through triggers in vector A
  for j in xrange(len(B)):
    # Iterate through triggers in vector B
    for i in xrange(len(A)):
      # Check if the difference is less than the time window
      if(np.abs(A[i] - B[j]) <=P):
	# If number of triggers is more than maxCoinc, exit.
	if(n<maxCoinc):
	  # Add those indexes to the list
	  C1.append(i)
	  C2.append(j)
	  n+=1
	else:
	  print "! Max num coincidences exceeded. Not recording further.\n"
	  cont=False
	  break
    if(cont==False):
      break
  C1 = np.asarray(C1)
  C2 = np.asarray(C2)
  return [C1, C2]

#def mfindcoinc(maxCoinc, A, B, P):
  ## Auxilary function used by mcoinc() to find co-incidences between time vectors
  ## A and B which are usually segments of the larger H and X channel trigger time
  ## co-ordniates
  ##
  ## Inputs: 
  ## maxCoinc - the maximum number of co-incidences
  ## A  - segment of the vector containing the time coordinates (seconds) triggers from a channel
  ## B  - same as above
  ## P - size of the coincidence window
  
  #n=0
  #cont=True
  #C1 = np.array([], dtype=int)
  #C2 = np.array([], dtype =int)
  ## Iterate through triggers in vector A
  #for j in xrange(len(B)):
    ## Iterate through triggers in vector B
      ## Check if the difference is less than the time window
    #idx = np.where((np.abs(A - B[j]) <=P))[0]
	## If number of triggers is more than maxCoinc, exit.
    #if(n<maxCoinc):
      ## Add those indexes to the list
      #C1 = np.append(C1, idx)
      #C2 = np.append(C2, j)
      #n+=1
    #else:
      #print "! Max num coincidences exceeded. Not recording further.\n"
      #cont=False
      #break
  #return [C1, C2]
        

def roundtopowertwo(segments, roundPower):
  numRound =  np.power(2.0, np.ceil(np.log2(segments)))
  
  numRoundMin = np.power(2.0, np.ceil(np.log2(roundPower)))
  
  if(numRound < numRoundMin):
    numRound = numRoundMin
  return numRound

def printvar(*varagarin):
  for i in range(0, len(varagarin), 2):
    print "%s %s\n" %(varagarin[i], varagarin[i+1])      

      
      
      
      
      
      
      
      
    




  
  
  