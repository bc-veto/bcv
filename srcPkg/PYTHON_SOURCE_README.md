Python Source for the Bilinear Coupling Veto
============================================
Readme file for code rewritten by Sudarshan Ghonge (<sudu.ghonge@gmail.com>) in Python based on code written by 
P Ajith
Bernard Hall
Nairwita Majumdar
Ben Yuan
and a few others.

To make things simple, the naming of variables, data flow and structure of functions has been kept the same.
All essential functions have been placed in a single module

bcv_runscript.py
================
BCV executable call script which takes in arguments in the same format as the 'omegaveto' file.

vetoanalysis.py
===============
Contains the function to which is called from the 'bcv_runscript' file

bcv.py
======
This module contains the following functions which serve the same functionality as the function used in the Matlab source.

readData
resample2
linearCouplingCoeff
bilinearCouplingCoeff
interpolatetransfn
calCrossCorr
highpass
loadframecache
readframedata
mcoinc
mfindcoinc
roundtopowertwo
printvar

trigstruct.py
============
File which contains a structure to store information about triggers.


