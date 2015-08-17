#!/usr/bin/python
import os;
import sys;
import math;
import string;

ifo = 'H1';
chunk = 'S6D';
outFileName = 'B-vetoAnalysis_metaDataSummary_%s_%s.mat' %(ifo, chunk);

rootDir = '/archive/home/detchar/public_html/S6/BCV/weekly';
#rootDir = '/archive/home/detchar/public_html/S6/BCV/S6C/';

# subdirectories 
subDirVecH1_S6A = [\
  'H1_DARMERR_931035296_931640096_WEEKLY_webpage',\
  'H1_DARMERR_931640096_932244896_WEEKLY_webpage',\
  'H1_DARMERR_932244896_932849696_WEEKLY_webpage',\
  'H1_DARMERR_932849696_933454496_WEEKLY_webpage',\
  'H1_DARMERR_933454496_934059296_WEEKLY_webpage',\
  'H1_DARMERR_934059296_934675200_WEEKLY_webpage',\
  'H1_DARMERR_934675200_935798487_WEEKLY_webpage'];

subDirVecH1_S6B = [\
  'H1_DARMERR_937800015_938822415_WEEKLY_webpage',\
  'H1_DARMERR_938822415_940032015_WEEKLY_webpage',\
  'H1_DARMERR_940032015_941241615_WEEKLY_webpage',\
  'H1_DARMERR_941241615_942451215_WEEKLY_webpage',\
  'H1_DARMERR_942451215_943660815_WEEKLY_webpage',\
  'H1_DARMERR_943660815_944870415_WEEKLY_webpage',\
  'H1_DARMERR_944870415_946080015_WEEKLY_webpage',\
  'H1_DARMERR_946080015_947260815_WEEKLY_webpage'];

subDirVecH1_S6C = [\
  'S6C_H1_DARMERR_949449543_950659287_inserted_webpage',\
  'S6C_H1_DARMERR_950659143_951868887_inserted_webpage',\
  'S6C_H1_DARMERR_951868743_953078487_inserted_webpage',\
  'S6C_H1_DARMERR_953078343_954288087_inserted_webpage',\
  'S6C_H1_DARMERR_954287943_955497687_inserted_webpage',\
  'S6C_H1_DARMERR_955497543_956707287_inserted_webpage',\
  'S6C_H1_DARMERR_956707143_957916887_inserted_webpage'];

subDirVecH1_S6D = [\
  'H1_DARMERR_959126415_959731215_WEEKLY_webpage',\
  'H1_DARMERR_959731215_960336015_WEEKLY_webpage',\
  'H1_DARMERR_960336015_960940815_WEEKLY_webpage',\
  'H1_DARMERR_960940815_961545615_WEEKLY_webpage',\
  'H1_DARMERR_961545615_962150415_WEEKLY_webpage',\
  'H1_DARMERR_962150415_962755215_WEEKLY_webpage',\
  'H1_DARMERR_962755215_963360015_WEEKLY_webpage',\
  'H1_DARMERR_963360015_963964815_WEEKLY_webpage',\
  'H1_DARMERR_963964815_964569615_WEEKLY_webpage',\
  'H1_DARMERR_964569615_965174415_WEEKLY_webpage',\
  'H1_DARMERR_965174415_965779215_WEEKLY_webpage',\
  'H1_DARMERR_965779215_966384015_WEEKLY_webpage',\
  'H1_DARMERR_966384015_966988815_WEEKLY_webpage',\
  'H1_DARMERR_966988815_967593615_WEEKLY_webpage',\
  'H1_DARMERR_967593615_968198415_WEEKLY_webpage',\
  'H1_DARMERR_968198415_968803215_WEEKLY_webpage',\
  'H1_DARMERR_968803215_969408015_WEEKLY_webpage',\
  'H1_DARMERR_969408015_970012815_WEEKLY_webpage',\
  'H1_DARMERR_970012815_970617615_WEEKLY_webpage',\
  'H1_DARMERR_970617615_971222415_WEEKLY_webpage',\
  'H1_DARMERR_971222415_971827215_WEEKLY_webpage'];
        
subDirVecL1_S6A = [\
  'L1_DARMERR_931035296_931640096_WEEKLY_webpage',\
  'L1_DARMERR_931640096_932244896_WEEKLY_webpage',\
  'L1_DARMERR_932244896_932849696_WEEKLY_webpage',\
  'L1_DARMERR_932849696_933454496_WEEKLY_webpage',\
  'L1_DARMERR_933454496_934059296_WEEKLY_webpage',\
  'L1_DARMERR_934059296_934675200_WEEKLY_webpage',\
  'L1_DARMERR_934675200_935798487_WEEKLY_webpage'];

subDirVecL1_S6B = [\
  'L1_DARMERR_937800015_938822415_WEEKLY_webpage',\
  'L1_DARMERR_938822415_940032015_WEEKLY_webpage',\
  'L1_DARMERR_940032015_940636815_WEEKLY_webpage',\
  'L1_DARMERR_940636815_941241615_WEEKLY_webpage',\
  'L1_DARMERR_941241615_941846415_WEEKLY_webpage',\
  'L1_DARMERR_941846415_942451215_WEEKLY_webpage',\
  'L1_DARMERR_942451215_943056015_WEEKLY_webpage',\
  'L1_DARMERR_943056015_943660815_WEEKLY_webpage',\
  'L1_DARMERR_943660815_944265615_WEEKLY_webpage',\
  'L1_DARMERR_944265615_944870415_WEEKLY_webpage',\
  'L1_DARMERR_944870415_945475215_WEEKLY_webpage',\
  'L1_DARMERR_945475215_946080015_WEEKLY_webpage',\
  'L1_DARMERR_946080015_946684815_WEEKLY_webpage',\
  'L1_DARMERR_946684815_947289615_WEEKLY_webpage'];

if ifo == 'H1': 
	if chunk == 'S6A':
		subDirVec = subDirVecH1_S6A;
	elif chunk == 'S6B':
		subDirVec = subDirVecH1_S6B;
	elif chunk == 'S6C':
		subDirVec = subDirVecH1_S6C;
	elif chunk == 'S6D':
		subDirVec = subDirVecH1_S6D;
elif ifo == 'L1': 
	if chunk == 'S6A':
		subDirVec = subDirVecL1_S6A;
	elif chunk == 'S6B':
		subDirVec = subDirVecL1_S6B;
	elif chunk == 'S6C':
		subDirVec = subDirVecL1_S6C;
	elif chunk == 'S6D':
		subDirVec = subDirVecL1_S6D;

#cmd = "nohup matlab -nodisplay -r convertsqlite2mat('%s','%s','%s') > convertsqlite2mat_%s_%s.log" \
#	%(ifo, dbFileLoc, outFileName, ifo, chunk);
for subDir in subDirVec: 

	dbFileLoc = '%s/%s/results' %(rootDir, subDir);	
	cmd = "matlab -nodisplay -r convertsqlite2mat('%s','%s','%s')" %(ifo, dbFileLoc, outFileName);
	cmd = cmd.replace(",","\,");
	cmd = cmd.replace("'","\\'");
	cmd = cmd.replace(")","\)");
	cmd = cmd.replace("(","\(");
	cmd = cmd.replace("/","\/");
	os.system(cmd);
