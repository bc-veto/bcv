# -*-shell-script-*-
#
# wsearch function definition
#
# Main Omega Pipeline function to search data for events
#
# usage: wsearch startTime stopTime
#
#   startTime           gps start time of analysis
#   stopTime            gps stop time of analysis
#
# variables:
#
#   parameterFile       parameter file
#   frameCacheFile      frame cache file
#   outputDirectory     directory to write results
#   logFile             alternate log file location
#   useTempDir          if set, use temporary staging directory
#   debugLevel          level of debug output
#
#   TMPDIR              location for temporary directories
#
# If the specified stop time is less than the specified start time, it is
# instead assumed to be the duration of the analysis.  Non-integer start and
# stop times are truncated to the nearest integer.
#
# If no parameter file is specified, wsearch looks for the file 'parameters.txt'
# in the current directory.  If that file is not found, or it is specified as
# '@online', and if the requested start time is less than 5 days old, wsearch
# selects the @online parameter set if it exists.  Other predefined parameter
# sets of likely interest are @H1, @H2, @L1, @H1H2, and @H1H2L1.
#
# Similarly, if no frame cache file is specified, wsearch looks for the file
# 'framecache.txt' in the current directory.  If that file is not found, it
# selects between @online and @default frame cache files based on the requested
# start time.  Other predefined frame cache files of likely interest are @S1,
# @S2, @S3, @A3, @S4, @A4, and @S5.  Alternatively, if the specified frame cache
# file argument points to a directory instead of a file, a temporary frame cache
# file will automatically be generated for the specified directory.
#
# The resulting triggers are collected into single files and stored in the
# following outputDirectory structure:
#
#   segments/
#     <startTime>-<stopTime>/
#       log.txt
#       livetime.txt
#       <channelName1>.txt
#       <channelName2>.txt
#       ...
#
# These files will contain redundant overlapping triggers as well as overlapping
# livetime segments.
#
# If an outputDirectory path is specified, then all results will be
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
# in the segment specific output directory.  This file contains the hostname of
# the machine where the job is running, and the path to the temporary segment
# specific output directory where the results are being written.  This file will
# be removed if the job exits cleanly.  If an existing lock file is found when a
# job starts, it will exit with an error.  Manual cleanup of the old temporary
# directory and lock file is then required.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>
# Jameson Graef Rollins <jrollins@phys.columbia.edu>

# $Id: wsearch.sh 1947 2009-07-06 17:45:53Z jrollins $

################################################################################
#                                                                              #
# cleanup                                                                      #
#                                                                              #
# Retrieves results from the temporary output directory, removes               #
# temporary and lock files, and exits.                                         #
#                                                                              #
################################################################################

cleanup() {

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

    # remove any temporary frame cache file
    if [ "${temporaryFrameCacheFile:-false}" = "true" ]; then
        rm -f "${frameCacheFile}"
    fi

    # remove lock file
    log "removing lock... "
    rm -f "${lockFile}"
    log "done.\n"
}

################################################################################
#                                                                              #
# wsearch                                                                      #
#                                                                              #
# Primary function to run wsearch binary executable                            #
#                                                                              #
################################################################################

wsearch() {

# ----------------------------------------
# process timing arguments
# ----------------------------------------

# check for valid number of arguments
if [ $# -ne 2 ] ; then
    failure "usage: wpipeline search [options] start stop"
fi

# extract start/stop times and truncate to integer
time0=$(printf '%0.0f' $1)
time1=$(printf '%0.0f' $2)

# calculate start/stop times
if (( time1 > time0 )) ; then
    startTime=$(( time0 ))
    stopTime=$(( time1 ))
elif (( time1 < time0 )) ; then
    startTime=$(( time0 ))
    stopTime=$(( time0 + time1 ))
else
    failure "ambiguous stop time"    
fi

# segment string (10 digit GPS)
segmentString=$(printf '%010.0f-%010.0f' $startTime $stopTime)

# round start time down to nearest integer GPS second
integerStartTime=$(printf '%0.0f' $startTime)

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
outputDirectory=${outputDirectory:-"segments/${segmentString}"}

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
    # make sure the TMPDIR exists
    mkdir -p "$TMPDIR"

    # make temporary directory
    temporaryDirectory=$(mktemp -d -t "wsearch_${segmentString}_XXXXXX")
    # canonicalize directory name
    temporaryDirectory=$(cd "$temporaryDirectory" && pwd)

    log "using tempdir: %s\n" "$temporaryDirectory"

    # copy any preexisting files to the temporary directory
    if find "${outputDirectory}" -mindepth 1 -maxdepth 1 | egrep -q '*' ; then
	log "copying files to temporary directory... "
	(cd "${outputDirectory}" && find . | grep -v '/.matlab' | \
	    cpio -pdu --quiet "${temporaryDirectory}")
	log "done.\n"
    fi
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
echo "$(hostname -f):${temporaryDirectory}" >"${lockFile}"

# ----------------------------------------
# run the executable
# ----------------------------------------

# set HOME to be the temporary directory, to avoid Matlab collisions
HOME="$temporaryDirectory"
export HOME

# call the wsearch executable
"${LIB}"/wsearch_bin \
    "$startTime" "$stopTime" \
    "$parameterFile" "$frameCacheFile" \
    "$temporaryDirectory" "$debugLevel" \
    | tee -a "$logFile" >> "$logOut"

}
