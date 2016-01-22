#!/bin/bash

# this script creates symbolic links to files in multiple folders
# and merges them into a single folder for use with .ini file
# that may only be able to accomodate a single trigger file location

#mkdir /home/bernard.hall/triggers

# ln -s /gds-h1/dmt/triggers/H-KW_TRIGGERS/H-KW_TRIGGERS-10780/* /home/bernard.hall/triggers/

SOURCE_ROOT=""
DEST_ROOT="/home/sudarshan.ghonge/"
DEST_FOLDER="triggers/"
DEST_DIR=$DEST_ROOT$DEST_FOLDER

echo "The destination directory is : "$DEST_DIR

ITERATOR=1
COUNTER=1

while [ $ITERATOR = 1 ]
do

S_CHECK=1

while [ $S_CHECK = 1 ]
do
echo "Enter the trigger directory tag"$COUNTER
read SOURCE_TAG
SOURCE_DIR=$SOURCE_ROOT$SOURCE_TAG
echo "The source directory is : "$SOURCE_DIR
echo "Is this correct?"
read S_RESPONSE

if [ $S_RESPONSE = "y" -o $S_RESPONSE = "Y" -o $S_RESPONSE = "yes" -o $S_RESPONSE = "Yes" ]
then
        S_CHECK=0
        #echo "ITERATOR=1"
else
        S_CHECK=1
        #echo "ITERATOR=0"
fi


done

#working part
#SOURCE_DIR="/gds-h1/dmt/triggers/H-KW_TRIGGERS/H-KW_TRIGGERS-10780"
echo "Working..."
DIR_LIST=("${DIR_LIST[@]}" $SOURCE_DIR)
echo "Done!"

echo "Enter another?"
read CHOICE

if [ $CHOICE = "y" -o $CHOICE = "Y" -o $CHOICE = "yes" -o $CHOICE = "Yes" ]
then
	ITERATOR=1
	#echo "ITERATOR=1"
        COUNTER=$(( COUNTER+1 ))
else
	ITERATOR=0
	#echo "ITERATOR=0"

fi   

#echo $CHOICE
#echo $ITERATOR

done
echo "${DIR_LIST[@]}"

for DIRECTORY in "${DIR_LIST[@]}"
do
	echo $DIRECTORY
	ln -s $DIRECTORY/* $DEST_DIR
done
echo "Exiting script!"
