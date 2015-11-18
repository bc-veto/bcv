import numpy as np
import matplotlib.pyplot as plt
import argparse
import bcv
import subprocess

parser = argparse.ArgumentParser(description = 'Plot the central times, durations and central frequencies of omicron triggers')

parser.add_argument('trighfile', type=str, help='Name of the h trigger file')
parser.add_argument('trigxfile', type=str, help='Name of the x trigger file')
parser.add_argument('startTime', type=int, help='Start times of the triggers')
parser.add_argument('endTime', type=int, help='End times of the triggers')

args = parser.parse_args()

trighfile = args.trighfile
trigxfile = args.trigxfile

trigHData = np.loadtxt(trighfile).reshape(-1, 9)
trigXData = np.loadtxt(trigxfile).reshape(-1, 9)
startTime = args.startTime
endTime = args.endTime

COINC_TIME_WINDOW = 0.5
segLength = 3600
uniqueArgument = 'nonunique'
maxNumCoinc = len(trigXData)*len(trigHData)
folder = subprocess.check_output('pwd').split('\n')[0]
channelXName = trigxfile.split('_', 3)[3].split('.')[0]
channelHName = trighfile.split('_', 3)[3].split('.')[0]

segIdx = np.intersect1d(np.where(trigXData[:, 2] > startTime)[0], np.where(trigXData[:,2] < endTime)[0])
trigXData = trigXData[segIdx]
segIdh = np.intersect1d(np.where(trigHData[:, 2] > startTime)[0], np.where(trigHData[:,2] < endTime)[0])
trigHData = trigHData[segIdh]

[coincTrigH, coincTrigX] = bcv.mcoinc(maxNumCoinc, trigHData[:, 2], trigXData[:, 2], COINC_TIME_WINDOW, segLength, uniqueArgument)

plt.figure()
plt.subplot(2,1,1)
ax = plt.gca()
ax.errorbar(trigXData[:, 2], trigXData[:, 3], xerr=[trigXData[:,2] - trigXData[:, 0], trigXData[:, 1] - trigXData[:, 2]], fmt='o')
plt.plot(trigXData[coincTrigX:,2], trigXData[coincTrigX:,3], '*')
plt.xlabel('Central Times')
plt.ylabel('Central Freqs')
ax.set_title('Triggers of %s'%(channelXName))

plt.subplot(2,1,2)
ax = plt.gca()
ax.errorbar(trigHData[:, 2], trigHData[:, 3], xerr=[trigHData[:,2] - trigHData[:, 0], trigHData[:, 1] - trigHData[:, 2]], fmt='o')
plt.plot(trigXData[coincTrigH:,2], trigXData[coincTrigH:,3], '*')
plt.xlabel('Central Times')
plt.ylabel('Central Freqs')
ax.set_title('Triggers of %s'%(channelHName))

plt.savefig(folder + 'Freq_time_%s.png'%(channelXName))
plt.close()




