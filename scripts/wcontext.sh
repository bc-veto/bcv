#!/bin/sh
#
# wcontext function definition
#
# Shell script to add timing, detector status, data quality, and
# detector log information to Omega Scans.
#
# usage: wcontext [options] gpsTime contextFile
#
# options: includeSections
#
# wcontext returns a list of all detector science mode and version 99
# data quality segments that occur within +/- 10 seconds of the
# requested gps time (rounded down to the nearest integer second).
#
# By default all information is included, but the optional
# includeSections argument can be used to display only a subset of the
# available information.  The includeSections argument should be a
# space separated list containing one or more of the keywords timing,
# segments, and logs.
#
# wcontext is called by wscan or wupdate and is not typically run
# directly by end users.  wcontext requires the wget, curl, or lynx
# utility as well as the LIGOTools segexpr utilities.  It is assumed
# that the location of these utilities are already included in the
# PATH environment variable as set by the calling function.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>

# $Id: wcontext.sh 2092 2009-08-08 22:56:06Z jrollins $

wcontext() {


# parse command line arguments
triggerTime="$1"
outDir="$2"


if [ -z "$includeSections" ]; then
    includeSections="timing segments logs"
fi

# context file name
contextFile="${outDir}/context.html"

# round trigger time down to nearest integer
integerTriggerTime=`echo ${triggerTime} | sed -e s'|\..*$||'`

# window to test around trigger time
windowDuration=10

# determine event directory
eventDirectory=`dirname ${contextFile}`

# list of interferometers to process
interferometers="G1 H1 H2 L1 V1"

# identify segments location and download utility
#if [ -d ${CACHE}/segments/S5 ]; then
#  get=cat
#  segmentsDirectory=${CACHE}/segments/S5
if [ -x "`which wget 2>/dev/null`" ]; then
  get="wget -q -O -"
  segmentsDirectory="http://ldas-cit.ligo.caltech.edu/segments/S5"
elif [ -x "`which curl 2>/dev/null`" ]; then
  get="curl"
  segmentsDirectory="http://ldas-cit.ligo.caltech.edu/segments/S5"
elif [ -x "`which lynx 2>/dev/null`" ]; then
  get="lynx -dump"
  segmentsDirectory="http://ldas-cit.ligo.caltech.edu/segments/S5"
else
  get=""
  segmentsDirectory=""
fi

# initialize context file
rm -f ${contextFile}
touch ${contextFile}

log "creating context file '${contextFile}'..."

# if timing information requested
if [ -n "`echo ${includeSections} | grep timing`" ]; then

  # begin timing section
cat <<EOF >> ${contextFile}
<a name="Timing"></a>
<h3>
<input id="input_Timing" type="checkbox" checked onclick="toggleVisible('Timing');" />
Event time
</h3>
<div class="section" id="div_Timing" style="display: block;">
<table>
<tbody>
EOF

  # date format string
  formatString="<tr><td><b>%Z:</b></td><td>%a</td>"
  formatString="${formatString}<td>%Y</td><td>%b</td><td>%d</td>"
  formatString="${formatString}<td>%H:%M:%S</td></tr>"

  # store time domain setting
  TZold=${TZ:-unset}

  # write GPS time
  echo "<tr><td><b>GPS:</b></td><td colspan=5 align=left>" \
       "${triggerTime}</td></tr>" >>${contextFile}

  # write UTC time
  TZ=UTC
  export TZ
  gpstime -g ${triggerTime} +"${formatString}" >>${contextFile}

  # write LHO time
  TZ="America/Los_Angeles"
  export TZ
  gpstime -g ${triggerTime} +"${formatString}" >>${contextFile}

  # write LLO time
  TZ="America/Chicago"
  export TZ
  gpstime -g ${triggerTime} +"${formatString}" >>${contextFile}

  # write VIRGO time
  TZ="Europe/Rome"
  export TZ
  gpstime -g ${triggerTime} +"${formatString}" >>${contextFile}

  # end timing section
cat <<EOF >> ${contextFile}
</tbody>
</table>
</div>

EOF

  # reset time domain setting
  if [ ${TZold} = "unset" ]; then
    unset TZ
  else
    TZ=${TZold}
    export TZ
  fi

# end test for timing information
fi

# if detector state and data quality information requested
if [ -n "`echo ${includeSections} | grep segments`" ]; then

  # if event time prior to S5
  if [ ${integerTriggerTime} -lt 815153408 ]; then

    # report unavailability of detector state and data quality information
    echo '<h3>Detector state and data quality information are not available<h3><br />' \
      >>${contextFile}

  # if event time is in S6 era
  elif [ ${integerTriggerTime} -gt 900000000 ]; then

    if [ -z "${categoryDefiner}" ]; then

	cat <<EOF >> ${contextFile}
<h3>Please use -q option with category definer file or ligolw dq tools for DQ info.</h3>
 
EOF
       
    else

cat <<EOF >> ${contextFile}
<a name="DQ"></a>
<h3>
<input id="input_DQ" type="checkbox" checked onclick="toggleVisible('DQ');" />
<a href=dq.txt>DQ information</a>
</h3>
<div class="section" id="div_Timing" style="display: block;">
<table>
<tbody>

EOF

        ligolw_dq_query --active $integerTriggerTime --database | ligolw_print --table segment_definer --column ifos --column name --delimiter ":" | sort -k1 > activeflags.tmp

        # Report whether each IFO is in Science Mode
        grep H1:DMT-SCIENCE activeflags.tmp > check.tmp
        if [ -s check.tmp ]; then
          echo '<p> H1 is in science mode' >> ${contextFile}
        else
          echo '<p><FONT COLOR="FF0000"> H1 is NOT in science mode </FONT>' >> ${contextFile}
        fi

	grep L1:DMT-SCIENCE activeflags.tmp > check.tmp
	if [ -s check.tmp ]; then
          echo '<p> L1 is in science mode' >> ${contextFile}
	else
          echo '<p><FONT COLOR="FF0000"> L1 is NOT in science mode </FONT>' >> ${contextFile}
	fi

	grep V1:ITF_SCIENCEMODE activeflags.tmp > check.tmp
	if [ -s check.tmp ]; then
          echo '<p> V1 is in science mode' >> ${contextFile}
	else
          echo '<p><FONT COLOR="FF0000"> V1 is NOT in science mode </FONT>' >> ${contextFile}
	fi
	rm check.tmp

	# erase local copy of category definer file if one exists
        if [ -e vetofile.xml ]; then
	    rm vetofile.xml
	fi    
	${get} "${categoryDefiner}" >> vetofile.xml
	export LIGOLW_VETO_FILE=vetofile.xml

	# get active, defined and inactive dq flags
	~omega/opt/glue/bin/ligolw_dq_active_cats $integerTriggerTime > "${outDir}/flagslist.tmp"

	# Add definition of category definer file and boilerplate
	cat <<EOF > "${outDir}/dq.txt"
# Data quality information obtained for GPS time ${integerTriggerTime} using category definer file:
# ${categoryDefiner}
# Information obtained $(date)
# Note that active flag definitions will vary depending on analysis.  The above flag states are valid only for the defined category definer file.
# <status> <category> <ifo>:<flag name>

# Active Data Quality Flags:
EOF
	grep '^+' "${outDir}/flagslist.tmp" | sort -k2 >> "${outDir}/dq.txt"

cat <<EOF >> ${outDir}/dq.txt

# Undefined Data Quality Flags:
EOF
        grep '^!' "${outDir}/flagslist.tmp" | sort -k2 >> "${outDir}/dq.txt"

cat <<EOF >> "${outDir}/dq.txt"

# Defined but Inactive Data Quality Flags:
EOF
        grep '^-' "${outDir}/flagslist.tmp" | sort -k2 >> "${outDir}/dq.txt"

    #End if for check on category definer file
    fi

cat <<EOF >> ${contextFile}
</tbody>
</table>
</div>

EOF

  # if segment information is not available
  elif [ -z "${segmentsDirectory}" ]; then

    # report unavailability of detector state and data quality information
    echo '<h3>Detector state and data quality information are not available<h3><br />' \
      >>${contextFile}

  # otherwise report detector state and data quality information
  else

    # download and parse segment information
    rm -f ${eventDirectory}/science_segments.txt
    rm -f ${eventDirectory}/dataquality_segments.txt
    LC_ALL=C
    export LC_ALL
    for interferometer in ${interferometers}; do
      ${get} ${segmentsDirectory}/${interferometer}/science_segments.txt \
        2>/dev/null | \
        awk ' !/^#/ { print $2, $3 } ' \
        >>${eventDirectory}/${interferometer}_science_segments.txt
      ${get} ${segmentsDirectory}/${interferometer}/injection_segments.txt \
        2>/dev/null | \
        awk ' !/^#/ { print $2, $3 } ' \
        >>${eventDirectory}/${interferometer}_science_segments.txt
      sort -g ${eventDirectory}/${interferometer}_science_segments.txt \
        >${eventDirectory}/${interferometer}_science_segments.txt.sorted
      mv ${eventDirectory}/${interferometer}_science_segments.txt.sorted \
         ${eventDirectory}/${interferometer}_science_segments.txt
      segexpr "union(${eventDirectory}/${interferometer}_science_segments.txt)" \
        2>/dev/null | \
        sed -e "s|^|${interferometer} Science_Mode |" \
        >>${eventDirectory}/science_segments.txt
      rm -f ${eventDirectory}/${interferometer}_science_segments.txt
      ${get} ${segmentsDirectory}/${interferometer}/dq_segments.txt \
        2>/dev/null | \
        awk ' $2 == 99 && $5 == 1 \
              { print "'${interferometer}'", $1, $3, $4 } \
              NR == 2 \
              { print "'${interferometer}' Not_Available", $2, "1e11" } ' \
        >>${eventDirectory}/dataquality_segments.txt
    done
    grep "Not_Available" ${eventDirectory}/dataquality_segments.txt \
      >>${eventDirectory}/science_segments.txt || true

    # determine time of last segment update
    updateTime=`awk ' /Not_Available/ { print $3 } ' \
                    ${eventDirectory}/science_segments.txt | \
                  sort -g | uniq | tail -n 1`
    if [ "$updateTime" ] ; then
	updateTimeString=`TZ=UTC gpstime -g ${updateTime} +"%a %Y %b %d %H:%M:%S %Z"`
    else
	updateTimeString='???'
    fi

    # begin detector state table
cat <<EOF >> ${contextFile}
<a name="Segments"></a>
<h3>
<input id="input_Segments" type="checkbox" checked onclick="toggleVisible('Segments');" />
Detector state
</h3>
(as of ${updateTimeString})
<div class="section" id="div_Segments" style="display: block;">
<table>
<tbody>
<tr>
<th>detector</th>
<th>state</th>
<th>start time</th>
<th>relative</th>
<th>stop time</th>
<th>relative</th>
</tr>
EOF

    # insert science mode segments
    cat ${eventDirectory}/science_segments.txt | \
      sed -e "s|$| ${integerTriggerTime}|" | \
      awk ' ($5 >= $3 - '${windowDuration}') && ($5 < $4 + '${windowDuration}') \
            { print $1, $2, $3, $3 - $5, $4, $4 - $5, $5 < $3 || $5 >= $4 } ' | \
      sort -k 7,7 -k 1,3 | \
      awk ' ($7 == 0) { printf "<tr>" } \
            ($7 == 1) { printf "<tr class=\"shade\">" } \
            { printf "<td>%s</td><td>%s</td><td>%d</td><td>%+d</td><td>%d</td><td>%+d</td></tr>\n", \
            $1, $2, $3, $4, $5, $6 } ' | \
      sed -e 's|Available</td>.*$|Available</td><td>-</td><td>-</td><td>-</td><td>-</td></tr>\n|' | \
      sed -e 's|_| |g' \
      >>${contextFile}

    # end detector state table
cat <<EOF >> ${contextFile}
</tbody>
</table>
</div>

EOF

    # begin data quality table
cat <<EOF >> ${contextFile}
<a name="Data_Quality"></a>
<h3>
<input id="input_Data_Quality" type="checkbox" checked onclick="toggleVisible('Data_Quality');" />
Data quality flags
</h3>
(as of ${updateTimeString})
<div class="section" id="div_Data_Quality" style="display: block;">
<table>
<tbody>
<tr>
<th>detector</th>
<th>flag</th>
<th>start time</th>
<th>relative</th>
<th>stop time</th>
<th>relative</th>
</tr>
EOF
    # insert data quality segments
    cat ${eventDirectory}/dataquality_segments.txt | \
      grep -v 'CBC_1YR_INJ_CHALLENGE' | \
      sed -e "s|$| ${integerTriggerTime}|" | \
      awk ' ($5 >= $3 - '${windowDuration}') && ($5 < $4 + '${windowDuration}') \
            { print $1, $2, $3, $3 - $5, $4, $4 - $5, $5 < $3 || $5 >= $4 } ' | \
      sort -k 7,7 -k 1,3 | \
      awk ' ($7 == 0) { printf "<tr>" } \
            ($7 == 1) { printf "<tr class=\"shade\">" } \
            { printf "<td>%s</td><td>%s</td><td>%d</td><td>%+d</td><td>%d</td><td>%+d</td></tr>\n", \
            $1, $2, $3, $4, $5, $6 } ' | \
      sed -e 's|Available</td>.*$|Available</td><td>-</td><td>-</td><td>-</td><td>-</td></tr>\n|' | \
      sed -e 's|_| |g' \
      >>${contextFile}

    # end data quality table
cat <<EOF >> ${contextFile}
</tbody>
</table>
</div>

EOF

    # remove temporary files
    rm -rf ${eventDirectory}/*_segments.txt

  # end test for availability of detector state and data quality information
  fi

# end test for detector state and data quality information
fi

# if detector log information requested
if [ -n "`echo ${includeSections} | grep logs`" ]; then

  # open detector log section
cat <<EOF >> ${contextFile}
<a name="Detector_Logs"></a>
<h3>
<input id="input_Detector_Logs" type="checkbox" checked onclick="toggleVisible('Detector_Logs');" />
Detector logs
</h3>
<div class="section" id="div_Detector_Logs" style="display: block;">
<p style="font-size: medium;">
EOF

  # store time domain setting
  TZold=${TZ:-unset}

  # links to LHO detector log
  TZ="America/Los_Angeles"
  export TZ
  LHOeventDay0=`gpstime -g ${triggerTime} +"%m/%d/%Y"`

  LHOeventDay1=`gpstime -g $(echo "$triggerTime + 86400" | bc -l) +"%m/%d/%Y"`
  LHOdetectorLog="http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?group=detector"
cat <<EOF >> ${contextFile}
LHO detector log:
<a href="${LHOdetectorLog}&date_to_view=${LHOeventDay0}">event day</a>,
<a href="${LHOdetectorLog}&date_to_view=${LHOeventDay1}">next day</a>
|
EOF

  # links to LLO detector log
  TZ="America/Chicago"
  export TZ
  LLOeventDay0=`gpstime -g ${triggerTime} +"%m/%d/%Y"`
  LLOeventDay1=`gpstime -g $(echo "$triggerTime + 86400" | bc -l) +"%m/%d/%Y"`
  LLOdetectorLog="http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?group=detector"
cat <<EOF >> ${contextFile}
LLO detector log:
<a href="${LLOdetectorLog}&date_to_view=${LLOeventDay0}">event day</a>,
<a href="${LLOdetectorLog}&date_to_view=${LLOeventDay1}">next day</a>
|
EOF

  # links to VIRGO detector log
  TZ="Europe/Rome"
  export TZ
  VIRGOeventDay0=`gpstime -g ${triggerTime} +"%m/%d/%Y"`
  VIRGOeventDay1=`gpstime -g $(echo "$triggerTime + 86400" | bc -l) +"%m/%d/%Y"`
  VIRGOeventDay2=`gpstime -g $(echo "$triggerTime + 172800" | bc -l) +"%m/%d/%Y"`
  VIRGOdetectorLog="https://pub3.ego-gw.it/logbook/index.php?area=logbook&ref=search"
cat <<EOF >> ${contextFile}
Virgo detector log:
<a href="${VIRGOdetectorLog}&datefrom=${VIRGOeventDay0}&dateto=${VIRGOeventDay1}">event day</a>,
<a href="${VIRGOdetectorLog}&datefrom=${VIRGOeventDay1}&dateto=${VIRGOeventDay2}">next day</a>
EOF

  # reset time domain setting
  if [ ${TZold} = "unset" ]; then
    unset TZ
  else
    TZ=${TZold}
    export TZ
  fi

  # close detector log section
cat <<EOF >> ${contextFile}
<p>
</div>
EOF

# end test for detector log information
fi

}
export -f wcontext
