# -*-shell-script-*-

# common variables and functions for Omega Pipeline interface

# Jameson Rollins <jrollins@phys.columbia.edu>

# $Id: common.sh 2108 2009-08-10 13:12:52Z jrollins $


### LOG
log() {
    if [ "$1" ] ; then
	printf "$@"
    else
	cat
    fi | sed 's/^/'"${LOG_PREFIX}"'/' >&2
}

### FAILURE
failure() {
    [ "$1" ] && echo "$1" >&2
    exit ${2:-'255'}
}

### SETUP_ENV
setup_env() {

    # setup the path
    PATH="$BIN":/usr/bin:/bin:"$PATH"
    export PATH

    # add the LAL library path
    LD_LIBRARY_PATH="$LAL"/lib:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH

    ### setup MATLAB standalone environment
    # mcr root path
    MCR_ROOT="${LIB}/mcr/${MCR_VNAME}"
    # mcr jre root path
    MCRJRE="${MCR_ROOT}"/sys/java/jre/${MATLAB_ARCH}/${JRE_VNAME}/lib/${JAVA_ARCH}

    LD_LIBRARY_PATH="${MCRJRE}":${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCRJRE}"/client:${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCRJRE}"/server:${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCRJRE}"/native_threads:${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCR_ROOT}"/sys/opengl/lib/${MATLAB_ARCH}:${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCR_ROOT}"/sys/os/${MATLAB_ARCH}:${LD_LIBRARY_PATH}
    LD_LIBRARY_PATH="${MCR_ROOT}"/runtime/${MATLAB_ARCH}:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH

    XAPPLRESDIR="${MCR_ROOT}"/X11/app-defaults
    export XAPPLRESDIR

    # variable to tell the matlab mcr to not lock the ctf archive
    MCR_INHIBIT_CTF_LOCK=1
    export MCR_INHIBIT_CTF_LOCK
}

### SETUP_TEMPDIR
# function to setup temporary directory
setup_tempdir() {
    local name="${1:-wpipeline}"

    # make sure the TMPDIR exists
    mkdir -p "$TMPDIR"

    # make temporary directory
    temporaryDirectory=$(mktemp -d -t "${name}_${eventString}_XXXXXX")
    # canonicalize directory name
    temporaryDirectory=$(cd "$temporaryDirectory" && pwd)

    log "using tempdir: %s\n" "$temporaryDirectory"
}

### XVFB_SETUP
# function to setup X virtual frame buffer (Xvfb) for functions
# that produce graphical output
xvfb_setup() {
    # variables
    local platform
    local lastTcpSocket
    local lastUnixSocket
    local lastDisplay
    local minimumDisplayNumber
    #displayNumber
    local xvfbTrials
    #xvfbProcessID
    local xvfbSuccess

    log "starting Xvfb... "

    platform=$(uname)

    # add path to Xvfb executable
    if [ -d /usr/X11R6/bin ]; then
	PATH=/usr/X11R6/bin:${PATH}
    elif [ -d /usr/openwin/bin ]; then
	PATH=/usr/openwin/bin:${PATH}
    fi
    export PATH

    # find last X11 tcp socket
    case ${platform} in
	'Linux')
	    lastTcpSocket=`netstat -lnt | grep tcp | \
		awk ' { print $4 } ' | grep ':60[0-9][0-9]$' | \
		sed -e 's|.*:60||' -e 's|^0||' | \
		sort -n | uniq | tail -1 || true`
	    ;;
	'SunOS')
	    lastTcpSocket=`netstat -aP tcp -f inet6 | grep '*.60' | \
		awk ' { print $1 } ' | grep '\.60[0-9][0-9]$' | \
		sed -e 's|.*\.60||' -e 's|^0||' | \
		sort -n | uniq | tail -1 || true`
	    ;;
	*)
	    failure "Unknown platform '$platform'."
	    ;;
    esac

    # find last X11 unix socket
    lastUnixSocket=`/bin/ls /tmp/.X11-unix/X* 2>/dev/null | \
        sed -e 's|.*X||' | sort -n | tail -1 || true`

    # find last X11 display number
    if [ "0${lastTcpSocket}" -gt "0${lastUnixSocket}" ]; then
	lastDisplay=${lastTcpSocket}
    else
	lastDisplay=${lastUnixSocket}
    fi
    
    # find unused X11 display number
    minimumDisplayNumber=50
    if [ -z "${lastDisplay}" ]; then
	displayNumber=${minimumDisplayNumber};
    elif [ ${lastDisplay} -lt ${minimumDisplayNumber} ]; then
	displayNumber=${minimumDisplayNumber};
    else
	displayNumber=`expr ${lastDisplay} + 1`
    fi

    # create virtual X display
    xvfbTrials=5
    while [ ${xvfbTrials} -gt 0 ]; do
	if [ "${platform}" = "Linux" ]; then
	    Xvfb :${displayNumber} -screen 0 1280x1024x24 -audit 0 -auth /dev/null -nolisten tcp 2>/dev/null &
	elif [ "${platform}" = "SunOS" ]; then
	    if [ `uname -r` = "5.8" ]; then
		Xvfb :${displayNumber} -audit 0 -auth /dev/null 2>/dev/null &
	    else
		Xsun :${displayNumber} -audit 0 -auth /dev/null +nkeyboard +nmouse \
                    -nolisten tcp -defdepth 24 -dev vfb 2>/dev/null &
	    fi
	fi
	xvfbProcessID=$!
	sleep 5
	xvfbSuccess=`/bin/ps -p ${xvfbProcessID} -o pid | sed -e 1d || true`
	if [ -n "${xvfbSuccess}" ]; then
	    break;
	fi
	xvfbTrials=`expr ${xvfbTrials} - 1`
	displayNumber=`expr ${displayNumber} + 1`
    done
    if [ -z "${xvfbSuccess}" ]; then
	failure "unable to open display."
    fi

    # point to virtual X display
    DISPLAY=":${displayNumber}.0"
    export DISPLAY

    # export variables
    export xvfbProcessID

    log "done.\n"
}

### XVFB_KILL
# kill the Xvfb process
xvfb_kill() {
    # the xvfb should be the only backgrounded job
    if jobs %1 &>/dev/null ; then
	log "stopping Xvfb... "
	kill %1
	wait
	log "done.\n"
    fi
    #kill ${xvfbProcessID} >/dev/null 2>&1;
}

### UTILITIES

### TIME2SCIENCE_RUN
# output the science run for a given GPS time
gps2science_run() {
    local time
    local run

    time="$1"

    if [ "$time" -lt 728955008 ]; then
	run="S1"
    elif [ "$time" -lt 751505200 ]; then
	run="S2"
    elif [ "$time" -lt 757711008 ]; then
	run="S3"
    elif [ "$time" -lt 792866944 ]; then
	run="A3"
    elif [ "$time" -lt 796144256 ]; then
	run="S4"
    elif [ "$time" -lt 815153408 ]; then
	run="A4"
    else
	run="S5"
    fi

    echo "$run"
}
