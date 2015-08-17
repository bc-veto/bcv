#!/bin/bash
#
# wpipeline: Interface to the Omega Pipeline burst search pipeline
#
# Copyright (C) 2008 Jameson Rollins <jrollins@phys.columbia.edu>
#
# This file is part of Omega Pipeline.
#
# Omega Pipeline is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# $Id$

################################################################################
#                                  variables                                   #
################################################################################

# fail on error
set -e

# set pipefail, so return code is that of failed command in pipe
set -o pipefail

# program name
PGRM=$(basename "$0" .sh)

# installation subdirectories
## DO NOT CHANGE ANY OF THESE ##
# these need to be set by various install targets
ROOT=__ROOT__
BIN="$ROOT"/__BIN__
export BIN
LIB="$ROOT"/__LIB__
export LIB
SHARE="$ROOT"/__SHARE__
export SHARE
CACHE="$ROOT"/__CACHE__
export CACHE
LAL=__LAL__
export LAL

# build info variables
VERSION=__VERSION__
MATLAB_VERSION=__MATLAB_VERSION__
MCR_VNAME=__MCR_VNAME__
MATLAB_ARCH=__MATLAB_ARCH__
JRE_VNAME=__JRE_VNAME__
JAVA_ARCH=__JAVA_ARCH__
FRAMEL_VERSION=__FRAMEL_VERSION__
LAL_VERSION=__LAL_VERSION__

# log output location
logOut=/dev/stderr

# make sure the useTempDir variable is not inherited from the environment
unset useTempDir

################################################################################
#                             function definitions                             #
################################################################################

usage() {
    cat <<EOF
Usage: $PGRM <subcommand> [options] [args]
Omega Pipeline command-line interface.

subcommands:
  search [options] <start> <stop>       search for triggers/events
    -p|--parameters <file>                parameter file
    -f|--framecache <file>                framecache file
  event [options] <center> <duration>   full analysis and plots of event time
    -p|--parameters <file>                parameter file
    -f|--framecache <file>                framecache file
    --event-file <file>                   event file to use for finding
                                          block times
  scan [options] <time>                 scan auxiliary channels for events
    -c|--configuration <file>             configuration file
    -f|--framecache <file>                framecache file
    -r|--report                           generate plots and html report
    -q|--veto-definer <file>              veto definer file for DQ reports
  configure-scan <configurationFile> <frameFile> [<frameFile> ...]
                                        create scan configuration file
  update-scan <scanDirectory>           update a scan directory
  properties                            produce property plots
    -p|--parameters <file>
    --triggers <file>
    --clusters <file>
    --segments <file>
    --injections <file>
    --livetime <file>
    --start <time>
    --stop <gpstime>
    --channel <name>
  post <function>                       exec post processing function
  online                                wonline subcommands
    setup <dir>                           setup a wonline directory
  version                               report version info
  help                                  this help

common options:
    -o|--outdir <dir>                   output directory
    -t[<dir>]|--tempdir[=<dir>]         temporary staging directory
    -l|--log <file>                     alternate log file location
    -s|--silent                         turn off logs to stderr
    -d|--debug <level>                  debug level

EOF
}

version() {
    printf "Omega Pipeline %s\n" \
	"$VERSION"
    printf "Matlab %s %s, MCR %s, JRE %s %s\n" \
	"$MATLAB_VERSION" \
	"$MATLAB_ARCH" \
	"$MCR_VNAME" \
	"$JRE_VNAME" \
	"$JRE_ARCH"
    printf "FrameL %s\n" \
	"$FRAMEL_VERSION"
    printf "LAL %s\n" \
	"$LAL_VERSION"
}

################################################################################
#                                  MAIN                                        #
################################################################################

### source common function definitions
source "$SHARE"/common.sh

if [ -z "$1" ] ; then
    failure "Type '$PGRM help' for usage."
fi

### parse command line
CMD="$1"
shift

TEMP=$(getopt --options p:f:c:rq:o:t::l:sd: --longoptions parameters:,framecache:,event-file:,configuration:,report,veto-definer:,triggers:,clusters:,segments:,injections:,livetime:,start:,stop:,channel:,outdir:,tempdir::,log:,silent,debug: -n "$PGRM" -- "$@")

if [ $? != 0 ] ; then
    usage
fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -p|--parameters)
            parameterFile="$2"
            shift 2
            ;;
        -f|--framecache)
            frameCacheFile="$2"
            shift 2
            ;;
	--event-file)
	    eventFile="$2"
	    shift 2
	    ;;
        -c|--configuration)
            configurationFile="$2"
            shift 2
            ;;
	-r|--report)
	    generateReport='true'
	    shift 1
	    ;;
	-q|--veto-definer)
	    categoryDefiner="$2"
	    shift 2
	    ;;
        --triggers)
            triggersFile="$2"
            shift 2
            ;;
        --clusters)
            clustersFile="$2"
            shift 2
            ;;
        --segments)
            segmentsFile="$2"
            shift 2
            ;;
        --injections)
            injectionsFile="$2"
            shift 2
            ;;
        --livetime)
            livetimeFile="$2"
            shift 2
            ;;
        --start)
            startTime="$2"
            shift 2
            ;;
        --stop)
            stopTime="$2"
            shift 2
            ;;
        --channel)
            channelName="$2"
            shift 2
            ;;
        -o|--outdir)
            outputDirectory="$2"
            shift 2
            ;;
        -t|--tempdir)
	    useTempDir="true"
	    if [ "$2" ] ; then
		export TMPDIR="$2"
	    fi
	    shift 2
	    ;;
	-l|--log)
	    logFile="$2"
	    shift 2
	    ;;
	-s|--silent)
	    logOut=/dev/null
	    shift 1
	    ;;
        -d|--debug)
            debugLevel="$2"
            shift 2
            ;;
        --)
            shift
            ;;
        *)
            break
            ;;
    esac
done

# prevent core dumps
ulimit -c 0

# explicity unset DISPLAY, since it seems to causes problems with the
# matlab binaries if it's set and not being captured with xvfb.  it is
# set explicitly by the xvfb function if needed.
unset DISPLAY

### source appropriate function definition and run
case $CMD in
    search)
        setup_env
        source "$SHARE"/wsearch.sh
        wsearch "$@"
        ;;
    event)
        setup_env
        source "$SHARE"/wevent.sh
        wevent "$@"
        ;;
    scan)
        setup_env
        source "$SHARE"/wscan.sh
        wscan "$@"
        ;;
    update-scan)
        source "$SHARE"/wupdate.sh
        wupdate "$@"
        ;;
    configure-scan)
        source "$SHARE"/wconfigure.sh
        wconfigure "$@"
        ;;
    properties)
        setup_env
        source "$SHARE"/wproperties.sh
        wproperties "$@"
        ;;
    online)
	SUB="$1"
	shift
	export WONLINE_PATH="${SHARE}/online"
	export PATH="$WONLINE_PATH":$PATH
	case $SUB in
	    setup)
		exec "$SHARE"/online/wonline_setup "$@"
		;;
	    supervisor)
		exec "$SHARE"/online/wonline_supervisor "$@"
		;;
	    search)
		exec "$SHARE"/online/wonline_search "$@"
		;;
	    *)
		failure "Unknown online subcommand '$SUB'.
Type '$PGRM help' for usage."
		;;
	esac
	;;
    post)
	function="$1"
	shift 1
	exec "$LIB"/postprocess/"$function" "$@"
	;;
    version)
	version
	;;
    'h'|'help'|'usage'|'-h'|'--help')
	usage
        ;;
    *)
	failure "Unknown command '$CMD'.
Type '$PGRM help' for usage."
        ;;
esac
