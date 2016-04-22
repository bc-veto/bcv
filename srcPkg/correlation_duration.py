import numpy as np
import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt
#import argparse

#parser = argparse.ArgumentParser(description = 'Script to generate correlation v/s duration plots of x triggers\n.
				 #The duration plotted against is of the x triggers')

#parser.add_argument('channelName', type=str, help='Name of the correlation statistic file which contains data such as correlation, duration, significance etc.')

#args = parser.parse_args()
#fileDir = args.channelName
if(len(sys.argv)<2):
  print 'Provide atleast one file name'

fileDir = sys.argv[1:]

channelName = fileDir[0].split('results/')[1].split('/')[0].split('+')[0].split('L1-')[1]
correlation = np.asarray([])
duration = np.asarray([])

for iFile in fileDir:
  corrData = np.loadtxt(iFile).reshape(-1, 15)
  correlation = np.append(correlation, corrData[:, 1])
  duration = np.append(duration, corrData[:, 9])
  
plt.figure()
plt.plot(duration, correlation, '.')
plt.title('Correlation vs trigger durations of %s channel' %(channelName))
plt.ylabel('correlation statistic $r$')
plt.xlabel('trigger durations')
plt.savefig('./corr_vs_duration_%s.png'%(channelName))
plt.close()
