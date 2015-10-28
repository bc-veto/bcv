import argparse
parser = argparse.ArgumentParser(description='HTML webpage creator for debug plots of the bilinear coupling veto pipeline')

parser.add_argument('results_directory', type=str, help='Results directory of the bcv runs')
parser.add_argument('webpage_directory', type=str, help='Directory where webpage is to be created')

args = parser.parse_args()

resultsDir = args.results_directory
webpageDir = args.webpage_directory
import glob
import dominate
from dominate.tags import *

# Store all debug plot folder paths for each pseudo channel
paths = glob.glob(resultsDir + '/*/*/debug_plots')

# Generate document title using the destination: The destination is usually the location of the webpage
doc_title=webpageDir.split('/')[-1] + ' debug_plots'

# Iterate through all pseudo channel
with dominate.document(title=doc_title) as doc:
  h1(doc_title)
  for ipath in paths:
    plot_info = ipath.split('/')
    # Get name of the pseudo channel
    pseudo_channel = plot_info[plot_info.index('results') + 1].split('-', 1)[1]
    h2(pseudo_channel)
    # Store the paths of all the timeshift folders in that particular pseudo channel
    timeShift_paths = glob.glob(ipath + '/*')
    # Iterate through each time shift
    for iTshift in timeShift_paths:
      tsfolder = iTshift.split('/')[-1]
      h3(tsfolder)
      # Store the paths of each of the plots in that time series folder
      plots = glob.glob(iTshift + '/*')
      # Iterate through each of the plots
      for iplot in plots:
	plotname = iplot.split('/')[-1]
	h4(plotname)
	h5('TimeSeries')
	div(img(src=iplot + 'TimeSeries.png'), _class='photo')
	h5('Specgram')
	div(img(src=iplot + 'Specgram.png'), _class='photo')
	

with open(webpageDir+ '/debug_plots.html', 'w') as f:
  f.write(doc.render())
    
    