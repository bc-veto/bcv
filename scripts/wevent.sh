# -*-shell-script-*-
#
# wevent function definition
#
# Omega Pipeline tool for analyzing an interesting event in detail.  wevent is a
# modified version of wsearch that graphically displays intermediate data
# products around a specified event time.
#
# usage: wevent eventTime duration
#
#   eventTime           gps center time of analysis
#   duration            plot duration
#
# variables:
#
#   parameterFile       parameter file
#   frameCacheFile      frame cache file
#   eventFile           event file in which to find event block info
#   outputDirectory     directory to write results
#   logFile             alternate log file location
#   useTempDir          if set, use temporary staging directory
#   debugLevel          level of debug output
#
#   TMPDIR              location for temporary directories
#
# If no parameter file is specified, wevent looks for the file 'parameters.txt'
# in the current directory.  If that file is not found, or it is specified as
# '@online', and if the requested start time is less than 5 days old, wevent
# selects the @online parameter set if it exists.  Other predefined parameter
# sets of likely interest are @H1, @H2, @L1, @H1H2, and @H1H2L1.
#
# Similarly, if no frame cache file is specified, wevent looks for the file
# 'framecache.txt' in the current directory.  If that file is not found, it
# selects between @online and @default frame cache files based on the requested
# start time.  Other predefined frame cache files of likely interest are @S1,
# @S2, @S3, @A3, @S4, @A4, and @S5.  Alternatively, if the specified frame cache
# file argument points to a directory instead of a file, a temporary frame cache
# file will automatically be generated for the specified directory.
#
# The results are placed in the following outputDirectory structure:
#
#   events/
#     <eventTime>
#       log.txt
#       index.html
#       summary.txt
#       summary.xml
#       ...
#
# If an outputDirectory path is specified, the all results will be
# placed there instead.
#
# If useTempDir is set, a temporary staging directory will be created
# in TMPDIR and all results will be written to this location during
# the analysis, and then transfered to the outputDirectory once the
# analysis is complete.  This is useful for running on clusters where
# the temporary directory can be specified on a local scratch disk on
# the cluster nodes.
#
# While the job is in progress, a lock file will be written to the file lock.txt
# in the event specific output directory.  This file contains the hostname of
# the machine where the job is running, and the path to the temporary event
# specific output directory where the results are being written.  This file will
# be removed if the job exits cleanly.  If an existing lock file is found when a
# job starts, it will exit with an error.  Manual cleanup of the old temporary
# directory and lock file is then required.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>
# Leo C. Stein <lstein@ligo.mit.edu>
# Jameson Graef Rollins <jrollins@phys.columbia.edu>

# $Id: wevent.sh 2313 2009-09-08 06:41:00Z jrollins $

################################################################################
#                                                                              #
# cleanup                                                                      #
#                                                                              #
# Retrieves results from the temporary output directory, removes               #
# temporary and lock files, and exits.                                         #
#                                                                              #
################################################################################

cleanup() {

    # kill Xvfb
    xvfb_kill

    # remove unneeded .matlab directory
    rm -rf "$temporaryDirectory"/.matlab

    # copy temporary files back to output directory and cleanup
    if [ "$useTempDir" -a -e "$temporaryDirectory" ] ; then
	log "copying files from temporary directory '%s'... " \
	    "$temporaryDirectory"
	(cd "${temporaryDirectory}" && find . | \
	    cpio -pdu --quiet "${outputDirectory}")
        rm -rf "${temporaryDirectory}"
	log "done.\n"
    fi

    # remove lock file
    log "removing lock... "
    rm -f "${lockFile}"
    log "done.\n"
}

################################################################################
#                                                                              #
# wevent                                                                       #
#                                                                              #
# Primary function to run wevent binary executable                             #
#                                                                              #
################################################################################

wevent() {

# ----------------------------------------
# process timing arguments
# ----------------------------------------

# check for valid number of arguments
if [ $# -ne 2 ] ; then
    failure "usage: wpipeline event [options] center duration"
fi

# extract event time arguments
eventTime="$1"
duration="$2"

# event string (10 digit GPS)
eventString=$(printf '%#020.9f' $eventTime)

# round event time down to nearest integer GPS second
integerStartTime=$(printf '%0.0f' $eventTime)

# start time of online data availability
onlineStartTime=$(( $(gpstime) - 5*24*3600 ))

# ----------------------------------------
# specify parameter file
# ----------------------------------------

# predefined parameter file directory
parameterDirectory="${SHARE}/parameters"

# set default parameter file
if [ -z "${parameterFile}" ]; then
    parameterFile=parameters.txt
    if [ ! -f "${parameterFile}" ]; then
        if [ "${integerStartTime}" -gt "${onlineStartTime}" ]; then
            parameterFile="@online"
        fi
    fi
fi

# handle predefined parameter files
if [ "$(echo ${parameterFile} | grep -q '^@')" ] ; then
    parameterFile=$(echo ${parameterFile} | \
        sed -e "s|^@\(.*\)$|${parameterDirectory}/\1.txt|")
fi

# ----------------------------------------
# specify frame cache file
# ----------------------------------------

# predefined frame cache file directory
frameCacheDirectory="${CACHE}/framecaches"

# set default framecache file
if [ -z "${frameCacheFile}" ]; then
    frameCacheFile=framecache.txt
    if [ ! -f "${frameCacheFile}" ]; then
        if [ "${integerStartTime}" -gt "${onlineStartTime}" ]; then
            frameCacheFile="@online"
        else
            frameCacheFile="@default"
        fi
    fi
fi

# frame cache file for online analysis
if [ "${frameCacheFile}" = "@online" ]; then
    if [ -d "/frames/full" ]; then
        frameCacheFile="/frames/full"
    elif [ -f "/virgoData/ffl/raw.ffl" ]; then
        frameCacheFile="/virgoData/ffl/raw.ffl"
    else
        frameCacheFile="@default"
    fi
fi

# frame cache file for default analysis
if [ "${frameCacheFile}" = "@default" ]; then
  frameCacheFile="@"$(gps2science_run "${integerStartTime}")
fi

# handle predefined frame cache files
if [ "$(echo ${frameCacheFile} | grep -q '^@')" ] ; then
    frameCacheFile=$(echo ${frameCacheFile} | \
        sed -e "s|^@\(.*\)$|${frameCacheDirectory}/\1.txt|")
fi

# generate frame cache file if requested
if [ -d "${frameCacheFile}" ]; then
    frameFileDirectory="${frameCacheFile}"
    dateString=$(date +%Y%m%d%H%M%S)
    frameCacheFile="/tmp/framecache_${dateString}_${integerStartTime}_$$.txt"
    rm -f "${frameCacheFile}"
    createframecache "${frameCacheFile}" "${frameFileDirectory}" >/dev/null
    temporaryFrameCacheFile=true
elif [ -n "$(echo ${frameCacheFile} | grep '\.ffl$')" ]; then
    fflFile="${frameCacheFile}"
    dateString=$(date +%Y%m%d%H%M%S)
    frameCacheFile="/tmp/framecache_${dateString}_${integerStartTime}_$$.txt"
    rm -f "${frameCacheFile}"
    convertfflcache "${fflFile}" "${frameCacheFile}"
    temporaryFrameCacheFile=true
fi

# ----------------------------------------
# define output directories
# ----------------------------------------

# set the default output directory
outputDirectory=${outputDirectory:-"events/${eventString}"}

# ----------------------------------------
# test lock
# ----------------------------------------

# set path to lock file
lockFile="${outputDirectory}/lock.txt"

# test for preexisting lock file
if [ -f "${lockFile}" ]; then
    failure "Lock file '${lockFile}' exists."
fi

# ----------------------------------------
# setup trap
# ----------------------------------------

# set trap to copy back output on exit
trap cleanup EXIT

# ----------------------------------------
# setup directories and log
# ----------------------------------------

# create output directory
mkdir -p "${outputDirectory}"

# canonicalize directory name
outputDirectory=$(cd "$outputDirectory" && pwd)

# create temporary directory if specified
if [ "$useTempDir" ] ; then
    setup_tempdir 'wevent'
else
    # else the "temporary" directory is the output directory
    temporaryDirectory="$outputDirectory"
fi

# set path to log file
logFile=${logFile:-"${temporaryDirectory}/log.txt"}

# ----------------------------------------
# create lock
# ----------------------------------------

# create lock file
printf "%s:%s\n" $(hostname -f) "${temporaryDirectory}" >"${lockFile}"

# ----------------------------------------
# add web support files
# ----------------------------------------

# copy html style sheet
cp "${SHARE}/misc/wstyle.css" "${outputDirectory}"

# copy icon
cp "${SHARE}/misc/logo/omega.logo.icon.png" "${outputDirectory}"

# note wevent is running in tempdir
if [ "$useTempDir" ] ; then
    cat <<EOF >"${outputDirectory}/index.html"
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>Event ${eventString}</title>
<link rel="icon" type="image/x-icon" href="omega.logo.icon.png" />
<link rel="stylesheet" type="text/css" href="wstyle.css" />
</head>
<body>
<h1 class="title"><a href="https://geco.phys.columbia.edu/omega" class="title">Omega Pipeline</a> Event</h1>
<h1>Event ${eventString}</h1>
<div class="main">
processing...<br />
$(date -u)<br />
$(hostname -f):${temporaryDirectory}<br />
</div>
</body>
</html>
EOF
fi

# ----------------------------------------
# setup display
# ----------------------------------------

# setup xvfb
xvfb_setup

# ----------------------------------------
# run the executable
# ----------------------------------------

# set HOME to be the temporary directory, to avoid Matlab collisions
HOME="$temporaryDirectory"
export HOME

# call the wevent executable
"${LIB}"/wevent_bin \
    "$eventTime" "$duration" \
    "$parameterFile" "$frameCacheFile" \
    "$temporaryDirectory" "$eventFile" \
    "$debugLevel" \
    | tee -a "$logFile" >> "$logOut"

}
