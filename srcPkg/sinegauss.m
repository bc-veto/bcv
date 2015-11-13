function [t, h] = sinegauss(dataLength, fs, t0, f0, phi0, Q, h_rss)
% 
% SINEGAUSS generate a sine gaussian waveform
% 
% usage: [t, h] = sinegauss(dataLength, fs, t0, f0, phi0, Q, h_rss)
% 
%     dataLength  - length of the data segment in secs
%     fs          - sampling frequency
%     t0          - time corresponding to the peak of the sine-gaussian (sec)
%     f0          - central frequency
%     phi0        - phase at t = t0
%     Q           - quality factor
%     h_rss       - RSS amplitude of the waveform
%
%     t           - time vector 
%     h           - waveform
%     
% P. Ajith, 29.10.2009
%
% $Id: mksinegauss.m,v 1.2 2006/04/06 14:04:43 ajith Exp $
   

% create a time vector
t  = 0:1/fs:dataLength-1/fs;

% generate the sine-gaussian 
tau       = Q/(sqrt(2)*pi*f0);
A         = h_rss*(2*f0^2/pi)^(1/4);
sinpart   = sin(2*pi*f0*(t-t0) + phi0);
gausspart = exp(-(t-t0).^2/(tau^2));
h         = A.*sinpart.*gausspart;

        
