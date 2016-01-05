import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys
import os
import bcv
parser = argparse.ArgumentParser(description = 'Plot the central frequency of x triggers vs h triggers')

parser.add_argument('trighfile', type=str, help='Name of the h trigger file')
parser.add_argument('trigxfile', type=str, help='Name of the x trigger file')

args = parser.parse_args()

trighfile = args.trighfile
trigxfile = args.trigxfile

# Force 2-D array shape
trigHData = np.loadtxt(trighfile).reshape(-1, 9)
trigXData = np.loadtxt(trigxfile).reshape(-1, 9)

# Time window of 0.5 s for co-incident triggers
COINC_TIME_WINDOW = 0.5
# Segment legnth in seconds to be fed to the mcoinc function
segLength = 3600
# Whether to have unique mapping b/w triggers
uniqueArgument = 'nonunique'
# Maximum number of coincidences
maxNumCoinc = len(trigXData)*len(trigHData)
# Get current working director
folder = os.popen('pwd').readlines()[0].split('\n')[0]

channelXName = trigxfile.split('_', 3)[3].split('.')[0]
channelHName = trighfile.split('_', 3)[3].split('.')[0]

[coincTrigH, coincTrigX] = bcv.mcoinc(maxNumCoinc, trigHData[:, 2], trigXData[:, 2], COINC_TIME_WINDOW, segLength, uniqueArgument)

if(len(coincTrigH)!=len(coincTrigX)):
  sys.exit('Inconsistent number of triggers')

plt.figure()
plt.plot(trigXData[:, 4][coincTrigX], trigHData[:, 4][coincTrigH], '.')
plt.xlabel('Central frequencies of X triggers')
plt.ylabel('Central frequencies of Y triggers')
plt.savefig(folder + 'Central_freqs_%s.png'%(channelXName))
plt.close()

