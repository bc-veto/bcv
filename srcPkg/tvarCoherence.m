function [ coherence ] = tvarCoherence( fchA,fchB,framePath,frameFile,start_time,end_time,increment_width,tslMin,tslMax, fs_multiplier, write_file )

%--------------------------------------
%if ~isnumeric(fs)
%	fs = str2num(fs);
%end
if ~isnumeric(start_time)
	start_time = str2num(start_time);
end
if ~isnumeric(end_time)
	end_time = str2num(end_time);
end
if ~isnumeric(increment_width)
	increment_width = str2num(increment_width);
end
if ~isnumeric(tslMin)
	tslMin = str2num(tslMin);
end
if ~isnumeric(tslMax)
	tslMax = str2num(tslMax);
end
if ~isnumeric(write_file)
	write_file = str2num(write_file);
end
if ~isnumeric(fs_multiplier)
        fs_multiplier = str2num(fs_multiplier);
end
%--------------------------------

duration = end_time - start_time;

[ts_1 fs1] = cohDataRead(framePath,frameFile,fchA,start_time,duration, write_file);

[ts_2 fs2] = cohDataRead(framePath,frameFile,fchB,start_time,duration, write_file);

% Pass fchA and fchB to tSlideCoherence to use created text files with
% write_file=1

% the lowest sampling rate will be chosen (higher will be downsampled in tSlideCoherence)
fs = 1;  % initialize sampling rate variable
rsflag = 1;

if fs1 > fs2
	fs = fs2;
	rsflag = 0;
end
if fs2 > fs1
	fs = fs1;
	rsflag = 0;
end

if rsflag
	fs = fs1;  % both sample rates were the same -> use first sample rate
end

% input sample rate will be reset to marked channel sample rates
%if (fs ~= fs1) && (rsflag == 1)  
%	disp('Requested sample rate is different from marked channel rate.  Resetting.')
%	fs = fs1;
%end

fs_string = strcat('Common sample rate:...',num2str(fs));
disp(fs_string)

coherence = tSlideCoherence( ts_1,ts_2,fs,tslMin,tslMax,increment_width,fs_multiplier,0 );

[m,n] = size(coherence);

for i = 1:n
	coherenceMat = cell2mat(coherence{i});
	dlmwrite(strcat('tests/coherence_',num2str(i),'.txt'),coherenceMat,' ');
end
