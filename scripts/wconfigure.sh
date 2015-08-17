# -*-shell-script-*-
#
# wconfigure function definition
#
# Shell script to automatically create a WSCAN configuration file for
# a given set of example frame file.
#
# usage: wconfigure configurationFile frameFile1 frameFile2 ...
#
#   configurationFile    path to output configuration file
#   frameFile1           path to first example frame file
#   frameFile2           path to second example frame file
#
# The FrChannels utility is used to extract the list of channel names
# and sample frequencies from the example frame file.  As a result,
# WCONFIGURE requires LIGOTOOLS.
#
# Each channel should only occur once in the set of all example frame
# files.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>

# $Id: wconfigure.sh 2072 2009-08-04 23:29:20Z bhughey $

wconfigure() {

# check for sufficiency command line arguments
if [ $# -lt 2 ]; then
    failure "wpipeline configure-scan configurationFile frameFile1 [frameFile2...]"
fi

# parse command line arguments
configurationFile="$1"
shift 1
frameFiles="$*"

# check for preexisting configuration file
if [ -f "${configurationFile}" ]; then
  failure "configuration file ${configurationFile} already exists"
fi

# initialize list of channels
channels=""

# begin loop over frame files
for frameFile in ${frameFiles}; do

  # check that frame file exists
  if [ ! -f ${frameFile} ]; then
    failure "frame file ${frameFile} does not exist"
  fi

  # extract frame file type from frame file name
  frameType=`basename ${frameFile} | sed -e 's|[^-]*-\([^-]*\)-.*|\1|'`

  # extract list of channel names and sample frequencies from frame file
  channels="${channels} `FrChannels ${frameFile} | \
            sed -e 's| |,|' -e "s|$|,${frameType}|"`"

# end loop over frame files
done

# filter out excluded detectors
excludedDetectors="HVE LVE"
for excludedDetector in ${excludedDetectors}; do
  channels=`echo ${channels} | sed -e 's| |\n|g' | \
            sed -e "/^${excludedDetector}/d"`
done

# sort by channel name
# unset die on failure so failed greps don't kill the script
set +e
sortedChannels=""
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                 sed -e 's| |\n|g' | \
                 sed -e 's|:| |' -e 's|,| |g' | \
                 sort -k 2,2 -k 1,1 | \
                 sed -e 's| |:|' -e 's| |,|g' | \
                 grep [^0]:LSC-STRAIN,`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 2,2 -k 1,1 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC-DARM_ERR,`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 2,2 -k 1,1 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC-DARM_CTRL,`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 2,2 -k 1,1 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC-AS_Q,`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 2,2 -k 1,1 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC-.*_EXC`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 2,2 -k 1,1 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC-.*_CAL | \
                grep -v [^0]:LSC-.*_EXC`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:LSC- | \
                grep -v [^0]:LSC-STRAIN, | \
                grep -v [^0]:LSC-DARM_ERR, | \
                grep -v [^0]:LSC-DARM_CTRL, | \
                grep -v [^0]:LSC-AS_Q, | \
                grep -v [^0]:LSC-.*_EXC | \
                grep -v [^0]:LSC-.*_CAL`"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:ASC- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:IOO- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:PSL- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:SUS- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:TCS- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:GDS- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:DAQ- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:SEI- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:OMC- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep [^0]:ISI- `"
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep -v 0: | \
                grep -v [^0]:LSC- | \
                grep -v [^0]:ASC- | \
                grep -v [^0]:IOO- | \
                grep -v [^0]:PSL- | \
                grep -v [^0]:SUS- | \
                grep -v [^0]:TCS- | \
                grep -v [^0]:GDS- | \
                grep -v [^0]:DAQ- | \
                grep -v [^0]:SEI- | \
                grep -v [^0]:OMC- | \
                grep -v [^0]:ISI- `"
                
sortedChannels="${sortedChannels} \
                `echo ${channels} | \
                sed -e 's| |\n|g' | \
                sed -e 's|:| |' -e 's|,| |g' | \
                sort -k 1,1 -k 2,2 | \
                sed -e 's| |:|' -e 's| |,|g' | \
                grep 0: `"
# reset die on failure
set -e

# initialize configuration file
touch "${configurationFile}"
cat <<EOF >>"${configurationFile}"
# Q Scan configuration file
# Automatically generated with wconfigure.sh
# by user ${USER} on $(date "+%F %T %Z")
# from sample frame files:
EOF

# add frame file comments
for frameFile in ${frameFiles}; do
  echo "#   ${frameFile}" >>"${configurationFile}"
done

# create default sections
cat <<EOF >>"${configurationFile}"

[Context,Context]

[Parameters,Parameter Estimation]

[Notes,Notes]

EOF

# initialize current subsystem name
currentSubsystemName=""

# begin loop over channels
for channel in ${sortedChannels}; do

  # display status
  [ "$channelName" ] && log "${channelName}\n"

  # extract channel name from channel list
  channelName=`echo ${channel} | sed -e 's|\([^,]*\),\([^,]*\),\([^,]*\)|\1|'`

  # extract sample frequency from channel list
  sampleFrequency=`echo ${channel} | sed -e 's|\([^,]*\),\([^,]*\),\([^,]*\)|\2|'`

  # extract frame file type from channel list
  frameType=`echo ${channel} | sed -e 's|\([^,]*\),\([^,]*\),\([^,]*\)|\3|'`

  # extract subsystem name from channel name
  case ${channelName} in
    *:LSC-STRAIN | *:LSC-DARM_ERR | *:LSC-DARM_CTRL | *:LSC-AS_Q)
      subsystemName="Gravitational"
      ;;
    *:LSC-*_EXC* | *:LSC-*_CAL*)
      subsystemName="Excitation"
      ;;
    *)
      subsystemName=`echo ${channelName} | sed -e 's|^\([^:]*:[^-]*\)-.*$|\1|'`
      ;;
  esac

  # if new subsystem
  if [ "x${subsystemName}x" != "x${currentSubsystemName}x" ]; then

    # update current subsystem name
    currentSubsystemName=${subsystemName}

    # extract detector name
    detectorName=`echo ${subsystemName} | sed -e 's|:.*||' \
                  -e 's|H0|LHO|' -e 's|L0|LLO|'`

    # do not scan vacuum equipment channels
    case ${detectorName} in
      HVE-* | LVE-*)
        continue
        ;;
    esac

    # determine section index and name from subsystem name
    case ${subsystemName} in
      Gravitational)
        sectionIndex=${subsystemName}
        sectionName="Gravitational wave data"
        ;;
      Excitation)
        sectionIndex=${subsystemName}
        sectionName="Calibration lines and injections"
        ;;
      *:LSC)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} length sensing and control"
        ;;
      *:ASC)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} alignment sensing and control"
        ;;
      *:IOO)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} input/output optics"
        ;;
      *:PSL)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} pre-stabilized laser"
        ;;
      *:SUS)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} suspension system"
        ;;
      *:TCS)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} thermal compensation system"
        ;;
      *:SEI)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} seismic isolation"
        ;;
      *:PEM)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} environment"
        ;;
      *:DAQ)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} data acquisition system"
        ;;
      *:GDS)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} global diagnostics system"
        ;;
      *:OMC)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} output mode cleaner"
        ;;
      *:ISI)
        sectionIndex=${subsystemName}
        sectionName="${detectorName} HAM6 ISI"
        ;;

      *:IFO | *:GPS | *:TID | *:TIM | *:FMC | *:GRB | *:DMT)
        # do not scan slow channels
        continue;
        ;;
      *)
        sectionIndex=${subsystemName}
        sectionName=${subsystemName}
        ;;
    esac

    # create new section
    cat <<EOF >>"${configurationFile}"
[${sectionIndex},${sectionName}]

EOF

  # end new subsystem
  fi

  # channel independent parameters
  searchQRange="[4 64]"
  searchMaximumEnergyLoss="0.2"
  plotFrequencyRange="[]"
  plotNormalizedEnergyRange="[0 25.5]"
  whiteNoiseFalseRate="1e-3"

  # begin setting parameters based on sample frequency
  case ${sampleFrequency} in

    16384)
      sampleFrequency="4096"
      searchTimeRange="64"
      searchWindowDuration="0.5"
      plotTimeRanges="[1 4 16]"
      ;;

    8192)
      sampleFrequency="4096"
      searchTimeRange="64"
      searchWindowDuration="0.5"
      plotTimeRanges="[1 4 16]"
      ;;

    4096)
      sampleFrequency="4096"
      searchTimeRange="64"
      searchWindowDuration="0.5"
      plotTimeRanges="[1 4 16]"
      ;;

    2048)
      sampleFrequency="2048"
      searchTimeRange="64"
      searchWindowDuration="0.5"
      plotTimeRanges="[1 4 16]"
      ;;

    1024)
      sampleFrequency="1024"
      searchTimeRange="64"
      searchWindowDuration="0.5"
      plotTimeRanges="[1 4 16]"
      ;;

    512)
      sampleFrequency="512"
      searchTimeRange="128"
      searchWindowDuration="1.0"
      plotTimeRanges="[1 4 16]"
      ;;

    256)
      sampleFrequency="256"
      searchTimeRange="256"
      searchWindowDuration="1.0"
      plotTimeRanges="[1 4 16]"
      ;;

    # do not scan slow channels
    *)
      continue
      ;;

  # end setting parameters based on sample frequency
  esac

  # set channel specific frequency ranges
  case ${channelName} in
    *:LSC-AS*_DAQ)
      sampleFrequency="64"
      searchTimeRange="1024"
      searchFrequencyRange="[0 Inf]"
      plotTimeRanges="[1 4 16]"
      ;;
    *:LSC-AS*_CORR_OUT_DAQ)
      sampleFrequency="64"
      searchTimeRange="1024"
      searchFrequencyRange="[0 Inf]"
      plotTimeRanges="[1 4 16]"
      ;;
    *:ASC-WFS*)
      sampleFrequency="64"
      searchTimeRange="1024"
      searchFrequencyRange="[0 Inf]"
      plotTimeRanges="[1 4 16]"
      ;;
    *:ASC-QPD*)
      sampleFrequency="64"
      searchTimeRange="1024"
      searchFrequencyRange="[0 Inf]"
      plotTimeRanges="[1 4 16]"
      ;;
    *:PEM-*_SEIS* | *:SEI-*_STS*)
      sampleFrequency="128"
      searchTimeRange="1024"
      searchFrequencyRange="[0 Inf]"
      plotTimeRanges="[8 64 512]"
      ;;
    *:SUS-*-OPLEV_*OUT)
      searchFrequencyRange="[0 Inf]"
      ;;
    *)
      searchFrequencyRange="[0 Inf]"
      ;;
  esac

  # set higher false rate for gravitational wave channels
  case ${channelName} in
    *:LSC-DARM_ERR)
      alwaysPlotFlag="1"
        ;;
    *:LSC-DARM_CTRL)
      alwaysPlotFlag="1"
        ;;
    *:LSC-AS_Q)
      alwaysPlotFlag="1"
        ;;
    *:LSC-STRAIN)
      alwaysPlotFlag="1"
        ;;
    *)
      alwaysPlotFlag="0"
        ;;
  esac

  # output channel configuation to file
  cat <<EOF >>"${configurationFile}"
{
  channelName:                 '${channelName}'
  frameType:                   '${frameType}'
  sampleFrequency:             ${sampleFrequency}
  searchTimeRange:             ${searchTimeRange}
  searchFrequencyRange:        ${searchFrequencyRange}
  searchQRange:                ${searchQRange}
  searchMaximumEnergyLoss:     ${searchMaximumEnergyLoss}
  whiteNoiseFalseRate:         ${whiteNoiseFalseRate}
  searchWindowDuration:        ${searchWindowDuration}
  plotTimeRanges:              ${plotTimeRanges}
  plotFrequencyRange:          ${plotFrequencyRange}
  plotNormalizedEnergyRange:   ${plotNormalizedEnergyRange}
  alwaysPlotFlag:              ${alwaysPlotFlag}
}

EOF

# end loop over channels
done

}
