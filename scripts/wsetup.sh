#!/bin/sh
#
# wsetup.sh
#
# Shell script to setup environment for Omega transform based tools.
#
# This script is sourced at the beginning of all Omega transform
# scripts and does not need to be sourced by the end user.
#
# This script sets up the PATH, LD_LIBRARY_PATH, MATLAB, MATLABPATH,
# and LIGOTOOLS environment variables to identify the location of
# MATLAB and LIGOTOOLS executable and libraries as well as other tools
# needed by Omega transform based applications.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>

# $Id: wsetup.sh 1171 2008-10-13 14:56:42Z jrollins $

# locations of required software
case `hostname` in
  ldas* | node* )
    MATLAB=/ldcg/matlab_r13
    LIGOTOOLS=/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  alterf | canopus )
    MATLAB=/apps/matlab
    LIGOTOOLS=/ligoapps/ligotools
    EXTRA_PATH=/ldcg/bin:/apps/workshop5/SUNWspro/bin
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=/home/shourov/matlab/curvefit/curvefit
    ;;
  fortress | padilla | prairie | snoqualmie | cresent | cougar | marble | evergreen )
    MATLAB=/export/apps/matlab6p5
    LIGOTOOLS=/export/apps4/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  decatur | basin | tuoro | delaronde | river | delta )
    MATLAB=/apps/stow/matlab-R13
    LIGOTOOLS=/apps/lscsoft/ligotools
    EXTRA_PATH=/home/shourov/olive/bin
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  defiance )
    MATLAB=/ldcg/matlab_r13
    LIGOTOOLS=/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  chestnut )
    MATLAB=/ldcg/matlab_r13
    LIGOTOOLS=/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  lancelot | enoki )
    MATLAB=/opt/matlab_6.5
    LIGOTOOLS=/net/emvogil-2/export/olive/packages/ligotools/ligotools
    EXTRA_PATH=/net/emvogil-2/export/olive/bin
    EXTRA_LD_LIBRARY_PATH=/net/emvogil-2/export/olive/packages/gcc/gcc-ldas-sparc-sun-solaris2.8/lib
    EXTRA_MATLABPATH=
    ;;
  *.psu.edu )
    MATLAB=/usr/global/lib
    LIGOTOOLS=/usr/global/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  *.uwm.edu )
    MATLAB=
    LIGOTOOLS=/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  *.virgo.infn.it | olmgr | olnode* )
    MATLAB=
    LIGOTOOLS=/users/shourov/home/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
  ccali* )
    MATLAB=/usr/local/matlab6p5
    LIGOTOOLS=/sps/virgo/USERS/shourov/ligotools
    EXTRA_PATH=
    EXTRA_LD_LIBRARY_PATH=
    EXTRA_MATLABPATH=
    ;;
esac

# determine platform name
platform=`uname`

# setup basic path
PATH=/bin:/usr/bin:/usr/local/bin
LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
if [ "${platform}" = "Linux" ]; then
  PATH=${PATH}:/usr/X11R6/bin
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/X11R6/lib
elif [ "${platform}" = "SunOS" ]; then
  if [ `uname -r` = "5.8" ]; then
    PATH=${PATH}:/usr/X11R6/bin:/usr/openwin/bin
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/X11R6/lib:/usr/openwin/lib
  else
    PATH=${PATH}:/usr/openwin/bin
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/openwin/lib
  fi
fi

# path to executable directory
executableDirectory=`dirname $0`

# add the executable directory to the executable path
PATH=${executableDirectory}:${PATH}

# add the executable directory to the matlab path
MATLABPATH=${executableDirectory}

# platform architecture name for Matlab libraries
if [ "${platform}" = "Linux" ]; then
  architecture=glnx86
elif [ "${platform}" = "SunOS" ]; then
  architecture=sol2
fi

# setup matlab r13 runtime libraries
if [ -d `dirname $0`/../lib/bin/${architecture} ]; then
  LD_LIBRARY_PATH=`dirname $0`/../lib/bin/${architecture}:${LD_LIBRARY_PATH}
fi

# setup matlab r13 environment
if [ -d ${MATLAB} ]; then
  PATH=${MATLAB}/bin:${PATH}
  LD_LIBRARY_PATH=${MATLAB}/bin/${architecture}:${LD_LIBRARY_PATH}
  LD_LIBRARY_PATH=${MATLAB}/extern/lib/${architecture}:${LD_LIBRARY_PATH}
  LD_LIBRARY_PATH=${MATLAB}/sys/os/${architecture}:${LD_LIBRARY_PATH}
  LD_LIBRARY_PATH=${MATLAB}/sys/opengl/lib/${architecture}:${LD_LIBRARY_PATH}
fi

# setup ligotools
if [ -n "${LIGOTOOLS}" ]; then
  PATH=${PATH}:${LIGOTOOLS}/bin
  MATLABPATH=${MATLABPATH}:${LIGOTOOLS}/matlab
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIGOTOOLS}/lib
fi

# include any extra path information
PATH=${PATH}${EXTRA_PATH:+:${EXTRA_PATH}}
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}${EXTRA_LD_LIBRARY_PATH:+:${EXTRA_LD_LIBRARY_PATH}}
MATLABPATH=${MATLABPATH}${EXTRA_MATLABPATH:+:${EXTRA_MATLABPATH}}

# prevent core dumps
ulimit -c 0

# export environment variables for use by subsequent commands
export PATH
export MATLAB
export LIGOTOOLS
export MATLABPATH
export LD_LIBRARY_PATH
