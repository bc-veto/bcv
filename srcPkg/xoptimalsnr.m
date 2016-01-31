function [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = xoptimalsnr(h,t0,fs,S,F0,dF,Fmin,Fmax)
%
% xoptimalsnr - Compute the SNR and characteristic frequency for a waveform
% in a specified noise background.  Simpler noise-independent measures of 
% the wave amplitude are also provided.
%
% This function is identical to OptimalSNR, except that the 
% plot option has been removed to make compiling easier.
%
%   [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = ... 
%       xoptimalsnr(h,t0,fs,S,F0,dF,Fmin,Fmax)
%
%   h      1 or 2-column array.  Waveform timeseries data.  Each column
%          holds the timeseries for one of the polarizations (plus or
%          cross).  The two polarizations are treated symmetrically, so it
%          does not matter which is the plus or the cross.
%   t0     Scalar.  Time at which the first waveform data point h(1,:) is
%          sampled.   
%   fs     Scalar.  Sampling rate (Hz) of waveform data.
%   S      Vector (optional).  Noise background one-sided POWER (not
%          amplitude) spectrum.  If supplied, it must cover at least the
%          range [Fmin, Fmax], and must be linearly sampled in ascending
%          order of frequency.  For noise-independent signal measures use
%          S=f0=df=[].  In this case the constant spectrum S(f)=2 will be
%          used (corresponding to a two-sided noise spectrum of unity).
%   F0     Scalar (optional).  Frequency at which S(1) is sampled.
%   dF     Scalar (optional).  Frequency spacing of the noise spectrum S.
%          THIS MUST EQUAL 1/(duration of timeseries h).
%   Fmin   Scalar (optional).  Minimum frequency (Hz) to include in 
%          frequency-domain calculations.
%   Fmax   Scalar (optional).  Maximum frequency (Hz) to include in 
%          frequency-domain calculations.
%
% Computations are done in the time domain (TD) and frequency
% domain (FD) using the energy distributions
%      p_TD = h(:,1).^2 + h(:,2).^2;
%      p_FD = 2(|FFT(h(:,1))|.^2 + |FFT(h(:,2))|.^2);  % for f>=0
% With these conventions, the output is
%   SNR    Signal-to-noise ratio of h in the given noise background, 
%          defined as 
%            SNR = (2 \int_Fmin^Fmax df p_FD./S).^0.5
%   h_rss  The root-sum-square amplitude (Hz^-0.5) of the waveform h:
%            h_rss = \int_Fmin^Fmax df p_FD
%   h_peak The maximum absolute amplitude of the waveform h:
%            h_peak = max(p_TD).^0.5
%   Fchar  Characteristic frequency (Hz) of h in the given noise  
%          background, defined as 
%                    \int_Fmin^Fmax df f p_FD./S
%            Fchar = ----------------------------------------
%                     \int_Fmin^Fmax df p_FD./S
%          where Fmin = max(f) and \tilde(h) is the FFT of h.  
%   bw     Effective bandwidth (Hz) of h in the given noise  
%          background, defined as 
%                 \int_Fmin^Fmax df (f-Fchar).^2 p_FD./S
%            bw = ---------------------------------------------------
%                  \int_Fmin^Fmax df p_FD./S
%          where Fmin = max(f) and \tilde(h) is the FFT of h.  
%   Tchar  Characteristic time at which the event occurs, defined as 
%                    \int dt t p_TD(t)
%            Tchar = -----------------
%                     \int dt p_TD(t)
%   dur    Duration of the signal, defined as 
%                    \int dt (t-Tchar).^2 p_TD(t)
%            Tchar = ----------------------------
%                     \int dt p_TD(t)
%
% Note that no windowing is used for computing FFTs; windowing and 
% averaging is not really sensible for a transient signal, since by
% definition the statistical properties of the signal are changing over
% time.  There may be problems for, e.g., band-passed noise bursts.
%
% Restriction: The waveform duration must equal 
% 1/(frequency resolution of f).
%
% Note that the SNR and h_rss measures are restricted to the frequency
% interval [Fmin,Fmax], h_peak is evaluated in the time domain (i.e., using
% data from the full frequency range. 
%
% Do "type FourierTransformConventions" for information on the conventions
% used for Fourier transforms.
% 
% $Id: xoptimalsnr.m 1941 2007-11-14 05:00:57Z jrollins $


%----- Preparations and checks.

%----- Optional arguments.
if ( (nargin<8) || isequal(Fmax,[]) )
    Fmax = floor(fs/2);
end
if ( (nargin<7) || isequal(Fmin,[]) )
    Fmin = 0;
end
%----- Is there an input noise spectrum?
if ( (nargin<6) || isequal(S,[]) || isequal(F0,[]) || isequal(dF,[]) )
    noise = 0;
else
    noise = 1;
end

%----- Define useful variables.
%----- Duration of tseries.
T = size(h,1)/fs;  
%----- Sample times.
t = t0+1/fs*[0:(size(h,1)-1)]';
%----- If no noise defined then make useful dummy noise vector of ones
%      covering [0,Nyquist] Hz. 
if (~noise)
    dF = 1/T;
    F0 = 0;
    S = 2*ones(floor(fs/2*T)+1,1);
end
%----- Vector of noise frequencies.
F = F0+[0:length(S)-1]'*dF;

%----- Error checks.
if (noise)
    if (dF*T ~= 1) 
        error('ERROR: Frequency resolution of noise spectrum does not match 1/duration of waveform.')
    end
    %----- Does vector of noise frequencies cover requested range?
    if ( (F(1)>Fmin) | (F(end)<Fmax)) 
        error('ERROR: Noise spectrum does cover desired frequency range.')
    end
end

%----- Convert everything to column vectors if needed.
if (size(S,2)>1)
    S = S';
end


%----- Frequency-domain calculations.

%----- Standard FFT.  Returns row/column vector according to input.
hf = 1/fs*fft(h);
%---------- Number of data points
N = size(hf,1);
%---------- Vector holding frequencies in usual screwy FFT order:
%           vector element:   [ 1  2  ...  N/2-1  N/2    N/2+1            N/2+2   ... N-1  N ]
%           frequency (df):   [ 0  1  ...  N/2-2  N/2-1  (N/2 or -N/2)   -N/2+1  ... -2   -1 ]
f = [ 0:N/2 , -N/2+1:-1 ]'/T;
%----- Take modulus, frequency components in [Fmin,Fmax]. 
%      THIS MEANS THAT WE ONLY KEEP POSITIVE FREQUENCIES.
%      INTEGRALS/SUMS OVER FREQUENCIES MAY NEED FACTOR 2 
%      TO ACCOUNT FOR NEGATIVE FREQUENCIES.
%---------- Signal:
k = find(f>=Fmin & f<=Fmax);
f = f(k);
p_FD = 2*sum(hf(k,:).*conj(hf(k,:)),2);
%---------- f=0 and f=Nyq bins should count for 1/2.
if (f(1)==0)
    p_FD(1) = 0.5*p_FD(1);
end
if (f(end)==floor(fs/2))
    p_FD(end) = 0.5*p_FD(end);
end    
%---------- Noise:
k = find(F>=Fmin & F<=Fmax);
F = F(k);
S = S(k);
%----- Error check
% if (isequal(f,F)==0)
%     [f(1) f(2) f(end)]
%     [F(1) F(2) F(end)]
%     error('Noise background not sampled at frequencies commensurate with time-series data.')
% end


%----- All set to do calculations.  Assured that f, F interchangable.

%----- SNR^2 versus frequency.
%SNRf = 4*(hm.^2)./S;
SNRf = 2*p_FD./S;

%----- SNR on [Fmin,Fmax].
SNR = (dF*sum(SNRf)).^0.5;

%----- Characteristic frequency.
Fchar = sum(f.*SNRf)./sum(SNRf);

%----- Characteristic bandwidth.  
bw = (sum((f-Fchar).^2.*SNRf)./sum(SNRf))^0.5;

% %------- Frequency interval containing central 50% of total |h(f)|^2.
% %------- Accumulation of |h(f)|^2 over frequency.
% chf2 = cumsum(SNRf)./sum(SNRf);
% %------- Select out unique points (necessary for linear interp)  
% [Z,I] = unique(chf2);
% chf2 = chf2(I);
% zf = f(I);
% f75 = interp1(chf2,zf,0.75);
% f25 = interp1(chf2,zf,0.25);
% bw = f75 - f25;

%----- RSS amplitude.
% h_rss = (2*dF*sum(hm.^2)).^0.5;
h_rss = (dF*sum(p_FD)).^0.5;


%----- TIME-DOMAIN calculations
p_TD = sum(h.^2,2);

%----- Peak amplitude.
% h_peak = max(abs(h));
h_peak = max(p_TD).^0.5;

%----- Characteristic time.
% Tchar = sum(t.*h.^2)./sum(h.^2);
Tchar = sum(t.*p_TD)./sum(p_TD);

%----- Duration.  
% dur = (sum((t-Tchar).^2.*h.^2)./sum(h.^2))^0.5;
dur = (sum((t-Tchar).^2.*p_TD)./sum(p_TD))^0.5;

% %------- Time containing 50% of total h(t)^2.
% %------- Accumulation of h(t)^2 over time.
% cht2 = cumsum(h.^2)/(h'*h);
% %------- Select out unique points (necessary for linear interp)  
% [Z,I] = UNIQUE(cht2);
% cht2 = cht2(I);
% zt = t(I);
% te = interp1(cht2,zt,0.75);
% ts = interp1(cht2,zt,0.25);
% dur = te - ts;


%----- Done
return
