function [x, t] = wsimulateddetectornoise(DET,T,fs,fmin,fmax,seed)
% wsimulateddetectornoise - Generate simulated Gaussian noise with spectrum 
% matching the design sensitivity curve for a specified gravitational-wave
% detector.
%  
%    [n, t] = SimulatedDetectorNoise(DET,T,fs,fmin,fmax)
%  
%    DET   String designating detector for which design sensitivity is 
%          desired.  Must be one of those recognized by the function SRD.
%    T     Scalar. Duration (s) of noise data stream.  T and fs should
%          satisfy T*fs=integer.
%    fs    Scalar. Sampling frequency (Hz) of data stream.  T and fs should
%          satisfy T*fs=integer.
%    fmin  Scalar. Minimum desired frequency (Hz).
%    fmax  Scalar. Maximum desired frequency (Hz).
%    seed  Optional scalar. Seed for randn generator used to produce the
%          noise data.  See randn.
%
%    n     Column vector of length round(fs*T).  Noise timeseries (strain).
%    t     Column vector of length round(fs*T).  Times (s) at which noise
%          vector is sampled, starting from t=0.
%  
% The design noise spectra are generated from the function SRD.  The noise
% power spectrum of the simulated data drops as f^(-2) outside the range 
% [fmin,fmax]. 
%    
% For information on the conventions used for discrete Fourier transforms,
% see FourierTransformConventions.m.
%
% -- original version: 
%      Patrick J. Sutton, 2005.04.09
%      psutton@ligo.caltech.edu
%  
%  $Id: SimulatedDetectorNoise.m 1013 2006-12-24 19:10:50Z psutton $

% ---- Checks.
error(nargchk(5,6,nargin));
if ~(T>0 && fs>0 && fmin>0 && fmax>0)
    error(['Input duration, sampling rate, and minimum and maximum ' ... 
        'frequencies must be positive.'])
end
if (fmax<=fmin)
    error('Maximum frequency must be > minimum frequency.')
end
if (fmax>fs/2)
    error('Maximum frequency must be <= Nyquist frequency fs/2.')
end
if (nargin==6)
    % if (~isscalar(seed))  % -- isscalar, isvector not in MatLab R13
    if (max(size(seed))>1)
        error('Seed value must be a scalar.')
    else
        % ---- Set state of randn. 
        randn('state',seed);
    end
end

% ---- Number of data points, frequency spacing.
N = round(fs*T);

% ---- Time at each sample, starting from zero.
t = [0:N-1]'/fs;

% ---- Make vector of positive frequencies up to Nyquist.
f = [1/T:1/T:fs/2]'; 

% ---- Get one-sided SRD (science requirements design) power spectrum (PSD).
PSD = SRD(DET,f);

% ---- Convert to two-sided power spectrum.
PSD = PSD/2;

% ---- Force noise spectrum to go to zero smoothly outside desired band.
k = find(f<fmin);
if (length(k)>0)
    PSD(k) = PSD(k(end)+1).*(f(k)/f(k(end))).^2;
end
k = find(f>fmax);
if (length(k)>0)
    PSD(k) = PSD(k(1)-1)*2./(1+(f(k)/f(k(1)-1)).^2);
end

% ---- Make white Gaussian random noise in frequency domain, at positive
%      frequencies (real and imaginary parts).
reXp = randn(length(f),1);
imXp = randn(length(f),1);

% ---- Color noise by desired amplitude noise spectrum.
Xp = ((T*PSD).^0.5).*(reXp + i*imXp)/2^0.5; 

% ---- Make noise at DC and negative frequencies, pack into vector in usual
%      screwy FFT order: 
%         vector element:   [ 1  2  ...  N/2-1  N/2    N/2+1            N/2+2   ... N-1  N ]
%         frequency (df):   [ 0  1  ...  N/2-2  N/2-1  (N/2 or -N/2)   -N/2+1  ... -2   -1 ]
%         F = [ 0:N/2 , -N/2+1:-1 ]'*df;
X = [ 0; Xp; conj(Xp(end-1:-1:1))];

% ---- Inverse fft back to time domain, casting off small imaginary part
%      (from roundoff error).
x = real(ifft(fs*X));

% ---- Done
return
