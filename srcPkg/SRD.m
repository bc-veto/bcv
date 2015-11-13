function [S, fstop] = SRD(IFO,f)
% SRD - Return science requirements design sensitivity curve for the
% specified gravitational-wave detector.
%
%  [S, fstop] = SRD(DET,f)
%
%  DET   String designating desired detector.
%        Must be one of 'AdvLIGO', 'GEO', 'LIGO', 'LIGO-LV', 'LHO', 'LLO',
%        'LISA', 'TAMA', 'VIRGO', 'Virgo-LV'.  Not case sensitive. 
%  f     Vector of positive frequencies (Hz) at which design sensitivity is
%        desired. 
%
%  S     Power spectral density (1/Hz) for requested detector and
%        frequencies.
%  fstop Minimum frequency (Hz) at which fit is "guaranteed" to be
%        accurate.  A warning message is printed if min(f)<fstop.
%
%  The fits used are from "Noise Power Spectrum of Initial
%  Interferometers", B. S. Sathyaprakash, Cardiff University, July 29 2000,
%  taken in turn from L. P. Grishchuk, V. M. Lipunov, K. A. Postnov, M. E.
%  Prokhorov, and B. S. Sathyaprakash, Physics-Uspekhi 71, 3 (2000).  The
%  'LIGO-LV' and 'VIRGO-LV' spectra are fits to the LIGO-Virgo Project Ib
%  simulated data sets.  'LIGO-LV' includes lines, while the others do not.
%  The 'LIGO-LV' reference frequency has been reset from 500 Hz to 300 Hz 
%  and other parameters tweaked to fit the simulated Project Ib Virgo data.
%  
%  -- original version:
%       Patrick J. Sutton, 2005.03.21
%       psutton@ligo.caltech.edu
%  
%  $Id: SRD.m 1952 2007-12-13 20:13:51Z jrollins $

% ---- Argument checks.
error(nargchk(2,2,nargin));
if ~ischar(IFO)
    error('Detector argument must be a string.')
end
% if ~isvector(f) % -- MatLab R13 does not have this function.
if (prod(size(f))~=length(f))
    error('Input frequencies must be a vector.')
end
if (min(f)<=0)
    error('Input frequencies must be positive.')
end

% ---- Recognize detector and compute noise spectrum.
switch upper(IFO)
    case 'ADVLIGO'
        fstop = 0;
        f0 = 215;
        x = f/f0;
        S = (1e-49)*( x.^(-4.14) -5*x.^(-2) + 111*(1-x.^2+x.^4/2)./(1+x.^2/2) );
    case 'GEO'
        fstop = 40;
        f0 = 150;
        x = f/f0;
        S = (1e-46)*( (3.4*x).^(-30) + 34*x.^(-1) + 20*(1-x.^2+x.^4/2)./(1+x.^2/2) );
    case {'LIGO', 'LHO', 'LLO'}
        fstop = 40;
        f0 = 150;
        x = f/f0;
        S = 9*(1e-46)*( (4.49*x).^(-56) + 0.16*x.^(-4.52) + 0.52 + 0.32*x.^2 );
    case 'LIGO-LV'
        %----- Design used by Chatterji in generating simulated data for 
        %      the LIGO-Virgo project Ia,b.  This design does not include
        %      the lines that were added to the simulated data.  fstop is 
        %      guessed.
        fstop = 40;
        seismicReference = 3.3e-21 * (30 ./ f).^14;
        thermalReference = 3.8e-23 * (100 ./ f).^2;
        shotReference = 1.13e-23 * sqrt(1 + (f / 90).^2);
        S = shotReference.^2 + thermalReference.^2 + seismicReference.^2;
    case 'LISA'
        fstop = 1e-5;
        f0 = 1e-3;
        x = f/f0;
        S = 420*(1e-46)*( (x/5.62).^(-14/3) + 1e3 + x.^2 );
    case 'TAMA'
        fstop = 75;
        f0 = 400;
        x = f/f0;
        S = 75*(1e-46)*( x.^(-5) + 13*x.^(-1) + 9 + 9*x.^2 );
    case 'VIRGO'
        fstop = 20;
        f0 = 500;
        x = f/f0;
        S = 3.24*(1e-46)*( (6.23*x).^(-5) + 2*x.^(-1) + 1 + x.^2 );
    case 'VIRGO-LV'
        fstop = 20;
        f0 = 500*0.6;  %-- reset from f0 = 500 to match LIGO-Virgo Project Ib simulated data
        x = f/f0;
        S = 3.24*(1e-46)*( (6.23*x).^(-5) + 2*x.^(-1) + 2 + 1.2*x.^2 );
    otherwise
        error(['Detector ' IFO ' not recognized.']);
end

% ---- Display warning message if asking for spectrum at frequencies below
%      the recommended cut-off.
if (min(f)<fstop)
    warning(...
        ['Requested frequency range extends below minimum frequency of ' ...
        num2str(fstop) ' Hz for SRD fit for ' IFO ' detector.']...
    )
end

% ---- Done.
return;

