# -*-shell-script-*-
#
# wupdate function definition
#
# Shell script to update timing, detector status, data quality, and
# detector log information of Omega Scans.
#
# usage: wupdate eventDirectory
#
# wupdate calls wcontext, which requires the wget, curl, or lynx
# utility.

# Shourov K. Chatterji <shourov@ligo.caltech.edu>

# $Id: wupdate.sh 1774 2009-04-27 06:20:20Z jrollins $

wupdate() {

source "$LIB"/wcontext.sh

eventDirectory="$1"

# determine index file
indexFile=${eventDirectory}/index.html

# validate index file
if [ ! -f ${indexFile} ]; then
    failure "Index file '${indexFile}' does not exist."
fi

# determine trigger time
triggerTime=`sed -e '/<title>/!d' \
                 -e 's|^[^0-9]*\([0-9]*.*[0-9]\)[^0-9]*$|\1|' \
               ${indexFile}`
if [ -z "${triggerTime}" ]; then
    failure "Could not determine event time."
fi

# determine context file
contextFile=${eventDirectory}/context.html

# initialize context file
rm -f ${contextFile}
touch ${contextFile}

# add timing and detector log information
wcontext ${triggerTime} ${contextFile}

# check for context markers
beginContextMarker=`grep '<!-- BEGIN CONTEXT -->' ${indexFile}`
if [ -z "${beginContextMarker}" ]; then
    failure "Could not locate BEGIN CONTEXT marker."
fi
endContextMarker=`grep '<!-- END CONTEXT -->' ${indexFile}`
if [ -z "${endContextMarker}" ]; then
    failure "Could not locate END CONTEXT marker."
fi

# insert context information
mv ${indexFile} ${indexFile}.tmp
sed -e '1,/<!-- BEGIN CONTEXT -->/!d' \
  ${indexFile}.tmp >>${indexFile}
echo "" >>${indexFile}
cat ${contextFile} >>${indexFile}
echo "" >>${indexFile}
sed -e '/<!-- END CONTEXT -->/,$!d' \
  ${indexFile}.tmp >>${indexFile}
rm -f ${indexFile}.tmp

}
