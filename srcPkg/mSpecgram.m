%
% MSPECGRAM mex file to compute spectrogram. 
%
% usage: [B,F,Tc] = mSpecgram(data,wind,nOverlap,nfft,fs)
% 
% data     : (vector) time-series data
% wind     : (vector) window function (of length nfft)
% nOverlap : overlaping samples between different FFTs
% nfft     : number of samples in one FFT
% fs       : sampling frequency of the data
%
% B        : complex FFT matrix (spectrogram) computed at positive
%            frequencies
% F        : frequencies at which the spectrogram is computed
% Tc       : center-times of the segments from which the spectrogram is computed
%            (note that the time Ts returned by the Matlab specgram function 
%            is the start time of the segment: Ts = Tc-nfft/2fs)
% 
% M Hewitson 15-02-06
%
% $Id: mSpecgram.m,v 1.2 2006/04/24 14:06:01 ajith Exp $
