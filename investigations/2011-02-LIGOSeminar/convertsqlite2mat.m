function convertsqlite2mat(ifo, dbFileLoc, outFileName)
% CONVERTSQLITE2MAT - convert the sqlite fines produced by the bilinear 
% veto code into a single mat file. 
% 
% usage: convertsqlite2mat(ifo, dbFileLoc, outFileName, iChunk)
% 
% P. Ajith, 4 April 2011 
% 
% $Id:$

% location and file name of the sqlite database file
if strcmp(ifo, 'H1') 
	dbFileVec_H1;
elseif strcmp(ifo, 'L1')
	dbFileVec_L1;
else
	error('ifo should be L1 or H1');
end

% if the mat file already exists, append the data 
if fopen(outFileName) >= 0
	load(outFileName)
	[nChunk, nChan] = size(MetaData);
	iChunk = nChunk+1;
else 
	iChunk = 1;
end

for iChan = 1:length(dbFileVec)
	dbFile = [dbFileLoc '/' dbFileVec{iChan}];
	if fopen(dbFile) >= 0 
		[MetaDataTmp, Data] = querydbfile(dbFile, []);
		File(iChunk,iChan).MetaDataTmp = MetaDataTmp;
		File(iChunk,iChan).MetaDataTmp.pseudoChannel = strrep(dbFileVec{iChan}, '.db', '');
		MetaData(iChunk,iChan) = File(iChunk,iChan).MetaDataTmp;
		clear File MetaDataTmp dbFile
	else 
		fprintf('... db file does not exist\n');
	end
end
    
% save the results into a mat file
save(outFileName, 'MetaData')
fprintf('... iChunk = %d. wrote mat file.\n', iChunk);
exit

