function [ coherence ] = tSlideCoherence( fchA,fchB,fs,start_time,end_time,increment_width, fs_multiplier, from_file)
%Performs coherence around a central time with an incrementally
%widening window for two input channels
%and returns an array of an array with a dataset from each
%successive window
%   An adaptation of code by Julius Orion Smith III,
%   from his book Mathematics of the DFT:
%   (http://www.dsprelated.com/dspbooks/mdft/Coherence_Function_Matlab.html)
%
% The sample frequency input to this function should be the lesser (if they
% are different) of the two input time series.  The time series with the
% greater sampling rate will be decimated to match that of the other
% series.

%-------------check parameters----------------------

if nargin == 3
    start_time = 0.1; %default beginning segment width
    end_time = 2.0;
    increment_width = 0.1;
    fs_multiplier = 1.0; 
    from_file = 0;
    
end
if nargin == 4
    end_time = 2.0; %default ending segment width
    increment_width = 0.1;
    fs_multiplier = 1.0;
    from_file = 0;

end
if nargin == 5
    increment_width = 0.1; %default increment width
    fs_multiplier = 1.0;
    from_file = 0;
end
if nargin == 6
    fs_multiplier = 1.0;
    from_file = 0;
end
if nargin == 7
    from_file = 0;
end

% -----------------------------------------------

if from_file  % expecting a single column of data (time series magnitude values)
        chA = load(fchA);
        chB = load(fchB);
else
        chA = fchA;
        chB = fchB;
end

if length(chA) > length (chB)
    %chA = downsample(chA,(length(chA)/length(chB)));
    chA = (decimate(chA,round(length(chA)/length(chB))))';
end

if length(chB) > length (chA)
    %chB = downsample(chB,(length(chB)/length(chA)));
    chB = (decimate(chB,round(length(chB)/length(chA))))';
end

N = length(chA); %the vectors should be the same length. Use channel A as a representative.  
central = N/2; %central time in samples

% -------------------------------------------------

i = 1; %coherence array index

for j = start_time:increment_width:end_time
    width = j; %full with of segment to be analyzed in seconds
    sample_width = (width * fs); %translation of width in seconds to width in samples
    lower = round(central - (0.5 * sample_width)); %lower bound in samples
    upper = round(central + (0.5 * sample_width)); %upper bound in samples

    if lower < 1   %do not allow bounds to go past the end of the array
        lower = 1;
    end

    if upper > (N - 1) %do not allow bounds to go past the end of the array
        upper = (N - 1);
    end
    x = chA(lower:upper,1); %select the proper length (in samples) of data for first channel
    y = chB(lower:upper,1); %select the proper length (in samples) of data for second channel
    nft = length(x);% nfft length.  Set to the length of the time series vectors. 
    
    if fs_multiplier == 0 % Set to zero, will find next power of 2 for nft
	p = nextpow2(nft);
	nft = 2^p;
    end

    %[cxy,f] = mscohere(x,y,[],[],128,fs);
    [cxy,f] = mscohere(x,y,[],[],nft,fs); % Do the work: nfft is set to the length of
    	% of time series; next power 2 higher (will be zero padded);

    vector_length = length(cxy);
    j_vec = zeros(vector_length,1);
    
    for k = 1:vector_length %create a vector filled with the value of j to act as the third column of 'coherence'
        j_vec(k) = j;
    end

    coherence{i} = {cxy,f,j_vec};
    i = i + 1; %(Post) increment the index
end
end 
