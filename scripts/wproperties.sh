#!/bin/sh
#
# wproperties function definition
#
# usage: wproperties parameterFile triggersFile clustersFile \
#                    segmentsFile injectionsFile livetimeFile \
#                    startTime stopTime channelName outputDirectory
#
# If not specified, the following default arguments are assumed.
#
#    parameterFile      ./parameters.txt
#    triggersFile       ./triggers.txt
#    clustersFile       ./clusters.txt
#    segmentsFile       ./segments.txt
#    injectionsFile     ./injections.txt
#    livetimeFile       ./livetime.txt
#    startTime          minimum time in livetime file
#    stopTime           maximum time in livetime file
#    channelName        channel name from parameter file
#    outputDirectory    ./properties
#    debugLevel         level of debug output

# Shourov K. Chatterji <shourov@ligo.caltech.edu>
# Jameson Graef Rollins <jrollins@phys.columbia.edu>

# $Id: wproperties.sh 1940 2009-07-06 03:46:51Z jrollins $

wproperties() {

# ----------------------------------------
# determine times
# ----------------------------------------

# determine start and stop times if not specified
if [ -z "${startTime}" ]; then
  startTime=`awk ' !/#/ { print $1 } ' ${segmentsFile} | sort -g | head -n 1`
fi
if [ -z "${stopTime}" ]; then
  stopTime=`awk ' !/#/ { print $2 } ' ${segmentsFile} | sort -g | tail -n 1`
fi

# determine channel name if not specified
if [ -z "${channelName}" ]; then
  channelName=`awk ' /channelNames/ ' ${parameterFile} | \
               sed -e "s|^[^']*'||" -e "s|'||g" -e "s| ||g" -e "s|}||g"`
fi

# time stamp
timeStamp=$(gpstime -g ${stopTime} +"%Y-%m-%d %H:%M:%S %Z")
duration=`expr ${stopTime} - ${startTime}`
hours=`expr ${duration} / 3600 | \
       awk ' { printf "%02d", $1 } '`
minutes=`expr \( ${duration} - ${hours} \* 3600 \) / 60 | \
         awk ' { printf "%02d", $1 } '`
seconds=`expr ${duration} - ${hours} \* 3600 - ${minutes} \* 60 | \
         awk ' { printf "%02d", $1 } '`
#timeStamp="${timeStamp} + ${hours}:${minutes}:${seconds}"
timeStamp="${timeStamp}"

# ----------------------------------------
# define output directories
# ----------------------------------------

# create output directory
mkdir -p "${outputDirectory}"

# set path to log file
logFile=${logFile:-"${outputDirectory}/log.txt"}

# ----------------------------------------
# make web pages
# ----------------------------------------

# copy html style sheet if one doesn't already exist
if [ ! -e "${outputDirectory}/wstyle.css" ] ; then
    cp "${SHARE}/wstyle.css" "${outputDirectory}"
fi

echo "$timeStamp" > "${outputDirectory}/timestamp.txt"

cat <<EOF >"${outputDirectory}/index.html"
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>Omega Online Properties: $startTime - $stopTime</title>
<link rel="stylesheet" href="wstyle.css" type="text/css">
</head>
<body>
<div class="content">
<h2>Omega Online Properties: $startTime - $stopTime</h2>
<div>
<h3>$timeStamp</h3><br>
<div>
<a href="unclustered_time_rate_5.png">
<img src="unclustered_time_rate_5.thumb.png" /></a>
<a href="clustered_time_rate_5.png">
<img src="clustered_time_rate_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_time_frequency_5.png">
<img src="unclustered_time_frequency_5.thumb.png" /></a>
<a href="clustered_time_frequency_5.png">
<img src="clustered_time_frequency_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_time_snr_5.png">
<img src="unclustered_time_snr_5.thumb.png" /></a>
<a href="clustered_time_snr_5.png">
<img src="clustered_time_snr_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_snr_rate_5.png">
<img src="unclustered_snr_rate_5.thumb.png" /></a>
<a href="clustered_snr_rate_5.png">
<img src="clustered_snr_rate_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_frequency_snr_5.png">
<img src="unclustered_frequency_snr_5.thumb.png" /></a>
<a href="clustered_frequency_snr_5.png">
<img src="clustered_frequency_snr_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_q_frequency_5.png">
<img src="unclustered_q_frequency_5.thumb.png" /></a>
<a href="clustered_q_frequency_5.png">
<img src="clustered_q_frequency_5.thumb.png" /></a>
</div>
<div>
<a href="unclustered_autocorrelogram_5.png">
<img src="unclustered_autocorrelogram_5.thumb.png" /></a>
<a href="clustered_autocorrelogram_5.png">
<img src="clustered_autocorrelogram_5.thumb.png" /></a>
</div>
</div>
</div>
</body>
</html>
EOF

# ----------------------------------------
# setup display
# ----------------------------------------

# setup xvfb
xvfb_setup

# ----------------------------------------
# run the executable
# ----------------------------------------

# call function to copy back output on exit
trap xvfb_kill EXIT

# call the wscan executable
"${LIB}"/wproperties_bin \
    "${parameterFile}" \
    "${triggersFile}" \
    "${clustersFile}" \
    "${segmentsFile}" \
    "${injectionsFile}" \
    "${livetimeFile}" \
    "${startTime}" "${stopTime}" "${channelName}" "${timeStamp}" \
    "${outputDirectory}" \
    "${debugLevel}" \
    | tee -a "${logFile}" > "$logOut"

}
