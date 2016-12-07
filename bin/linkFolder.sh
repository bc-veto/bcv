#!/bin/bash

# this script creates symbolic links to files in multiple folders
# and merges them into a single folder for use with .ini file
# that may only be able to accomodate a single trigger file location

#mkdir /home/bernard.hall/triggers

# ln -s /gds-h1/dmt/triggers/H-KW_TRIGGERS/H-KW_TRIGGERS-10780/* /home/bernard.hall/triggers/
DIRECTORY=$1
DEST_FOLDER=$2
echo "Source directory is " $DIRECTORY
echo "Destination directory is " $DEST_FOLDER
ln -s $DIRECTORY/* $DEST_FOLDER


