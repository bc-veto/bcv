function [t,hp,hc] = xmakewaveform(type,parameters,T,T0,fs,varargin)
% xmakewaveform - Create timeseries data containing a specified type of 
% simulated gravitational-wave signal.
%
% Note: This function is largely identical to MakeWaveform, except that it
% calls xoptimalsnr instead of OptimalSNR; the former has graphics commands 
% removed to make compiling easier.
%
%   [t,hp,hc] = xmakewaveform(type,params,T,T0,fs,'PropertyName', ...
%       'PropertyValue',...)
%
%   type    A char or cell array containing one of 
%            'CG'  cosine-Gaussians (special case of 'chirplet'), 
%      'chirplet'  Chirping sine/cosine-Gaussian,
%          'cusp'  waveform from a cosmic string cusp (see cuspgw.m),
%           'DFM'  Dimmelmeier-Font-Mueller supernova waveforms,
%          'DFMc'  Conditioned Dimmelmeier-Font-Mueller supernova 
%                  waveforms.  Same as DFM except detrended.
%            'DS'  damped sinusoids, 
%          'DS2P'  2-polarization damped sinusoids with smooth turn-on, 
%             'G'  Gaussians (special case of 'chirplet'), 
%       'Lazarus'  Lazarus project black-hole merger waveforms (analytic 
%                  approximation from Baker, Campanelli, Lousto, and 
%                  Takahashi, PRD 65 124012 2002. ),
%            'OB'  Ott-Burrows supernova waveforms.
%  'onecyclesine'  single cycle (phase -pi->0->pi) of a sine wave,
%            'SG'  sine-Gaussians (special case of 'chirplet'), 
%           'WNB'  Windowed, band-limited white-noise burst (2 independent 
%                  polarizations). 
%          'zero'  null signal (hp=0=hc),
%            'ZM'  Zwerger-Mueller supernova waveforms.
%   params  Array (double or cell), or tilde-delimited string (see below) 
%           containing parameters for the specified waveform type:
%             CG:  Double array [h_peak,tau,f0], where 'h_peak' is the peak   
%                  amplitude, 'tau' is the duration, and 'f0' is the  
%                  central frequency.
%                      hp = h_peak*cos(2*pi*(t-T0)*f0).*exp(-(t-T0).^2./tau.^2)
%                      hc = 0;
%      'chirplet'  Double array [hrss,tau,f0,alpha,delta], where hrss is
%                  the quadrature sum of the RSS amplitudes of the plus and
%                  cross polarizations, tau is the duration, f0 is the
%                  central frequency, alpha is the chirp parameter, and
%                  delta is the phase at the peak of the envelope.
%                     hp = real(h);
%                     hc = imag(h);
%                  where
%                     h = hrss*exp(...
%                             - (1-i*alpha)*(t-T0).^2./(4*tau.^2) ...
%                             + i*2*pi*(t-T0)*f0 + i*delta  ...
%                         )./(2*pi*tau^2).^(1/4);
%                  Note that alpha>0 gives a signal that chirps upward in
%                  frequency.  The signal time-frequency volume is
%                     V = (1+alpha^2)^0.5;
%                  and the bandwidth is
%                     bw = V/(4*pi*tau) = (1+alpha^2)^0.5/(4*pi*tau);
%           cusp:  Scalar f0, which is the cutoff frequency in Hz.
%            DFM:  Cell array [distance,designation], where 'distance' is 
%                  the distance to the supernova in kpc and 'designation' 
%                  is the waveform type.  Known waveform types are 
%                  'A1B1G1', 'A1B2G1', 'A1B3G1', 'A1B3G2', 'A1B3G3',
%                  'A1B3G5', 'A2B4G1', 'A3B1G1', 'A3B2G1', 'A3B2G2',
%                  'A3B2G4', 'A3B2G4_soft', 'A3B3G1', 'A3B3G2', 'A3B3G3',
%                  'A3B3G5', 'A3B4G2', 'A3B5G4', 'A4B1G1', 'A4B1G2',
%                  'A4B2G2', 'A4B2G3', 'A4B4G4', 'A4B4G5', 'A4B5G4',
%                  'A4B5G5'.
%           DFMc:  Same as DFM.  
%             DS:  Double array [h_peak,tau,f0], where 'h_peak' is the peak   
%                  amplitude, 'tau' is the duration, and 'f0' is the  
%                  central frequency:
%                         hp = h_peak*cos(2*pi*(t-T0)*f0).*exp(-(t-T0)/tau);  t>=T0
%                            = 0;  t<t0
%                         hc = h_peak*sin(2*pi*(t-T0)*f0).*exp(-(t-T0)/tau);  t>=T0
%                            = 0;  t<t0
%           DS2P:  Double array [h_peak,tau,f0,delta,ciota], where 'h_peak'    
%                  is the peak amplitude, 'tau' is the duration, 'f0' is  
%                  the central frequency, 'delta' is a phase, and 'ciota' 
%                  is the cosine of the inclination angle. 
%                         hp = h_peak*0.5*(1+ciota^2)*cos(2*pi*(t-T0)*f0+delta).*exp(-(t-T0)/tau);     t>=T0  
%                            = h_peak*0.5*(1+ciota^2)*cos(2*pi*(t-T0)*f0+delta).*exp((t-T0)/(tau/10)); t<T0  
%                         hc = h_peak*ciota*sin(2*pi*(t-T0)*f0+delta).*exp(-(t-T0)/tau);     t>=T0  
%                            = h_peak*ciota*sin(2*pi*(t-T0)*f0+delta).*exp((t-T0)/(tau/10)); t<T0  
%              G:  Double array [h_peak,tau], where 'h_peak' is the peak amplitude of  
%                  the waveform, 'tau' is the duration, following the 
%                  conventions of the S1 Bursts paper (gr-qc/0312056).
%                    hp = h_peak*exp(-(t-T0).^2/tau.^2);
%                    hc = 0;
%        Lazarus:  [mass,distance,ciota], where 'mass' is the binary system
%                  mass in solar masses, 'distance' is the distance to the
%                  system in Mpc, and 'ciota' is the cosine of the 
%                  inclination angle.  The inclination angle is the angle
%                  between the line of sight to the system and the system
%                  rotational axis.  The amplitudes of the two
%                  polarizations vary with ciota as 
%                    hp  \propto  1/2*(1+ciota^2);
%                    hc  \propto  ciota;
%             OB:  Cell array [distance,designation], where 'distance' is 
%                  the distance to the supernova in kpc and 'designation' 
%                  is the waveform type.  Known waveform types are 
%                  'e15', 'e20', 'm15b4', 'm20b4', 'm25b4V',
%                  's11A1000B0.1', 's11A1000B0.2', 's11A1000B0.3',
%                  's11A1000B0.4', 's11A1000B0.5', 's11A1000B0.6',
%                  's11A1000B0.7', 's11A1000B0.8', 's11A50000B0.1',
%                  's11A50000B0.2', 's11A50000B0.25', 's11A50000B0.3',
%                  's11A50000B0.4', 's11A50000B0.5', 's11A50000B0.6',
%                  's11A50000B0.7', 's11A500B0.1', 's11A500B0.2',
%                  's11A500B0.25', 's11A500B0.3', 's11A500B0.4',
%                  's11A500B0.5', 's11nonrot', 's15A1000B0.1',
%                  's15A1000B0.2', 's15A1000B0.3', 's15A1000B0.4',
%                  's15A1000B0.5', 's15A1000B0.6', 's15A1000B0.7',
%                  's15A1000B0.8', 's15A1000B0.9', 's15A1000B1.0',
%                  's15A50000B0.1', 's15A50000B0.2', 's15A50000B0.5',
%                  's15A50000B1.0', 's15A500B0.1', 's15A500B0.2',
%                  's15A500B0.25', 's15A500B0.3', 's15A500B0.4',
%                  's15A500B0.5', 's15A500B0.6', 's15A500B0.9',
%                  's15A500B1.0', 's15nonrot'.
%   onecyclesine:  Double array [h_peak,f0], where 'h_peak' is the peak 
%                  amplitude of the waveform, and 'f0' is the frequency of 
%                  the sine wave:
%                      hp = h_peak*sin(2*pi*(t-T0)*f0)  for  abs(t-T0) < 1/(2f0);  
%                      hp = 0  for  abs(t-T0) >= 1/(2f0);  
%                      hc = 0;
%             SG:  Double array [h_peak,tau,f0], where 'h_peak' is the peak   
%                  amplitude, 'tau' is the duration, and 'f0' is the  
%                  central frequency, following the conventions of the S1
%                  bursts paper (gr-qc/0312056).
%                      hp = h_peak*sin(2*pi*(t-T0)*f0).*exp(-(t-T0).^2./tau.^2); 
%                      hc = 0;
%            WNB:  Double array [hrss, fc, df, dt] where 'hrss' is the 
%                  root-sum-square signal amplitude, 'fc' is the desired
%                  central frequency, 'df' is the desired bandwidth, and
%                  'dt' is the width (duration) of the Gaussian envelope. 
%                  Each polarization is an independent sample of noise
%                  that is white over the band [max(0,fc-df),min(fc+df,fs)] 
%                  and zero elsewhere, modulated by the Gaussian envelope
%                      env = exp(-(t-T0).^2/2/dt.^2) 
%                  Note that certain combinations of parameters are not
%                  realizable.  For example, all signals must have
%                  df*dt>~0.1.  Note also that these waveforms are
%                  generated using random numbers, so a particular
%                  instantiation may not have precisely the desired
%                  duration or frequency content (usually good to a few x
%                  10%).  Use OptimalSNR to verify that the output waveform
%                  has the desired time-frequency properties.  
%           zero:  [] (empty).
%             ZM:  Cell array [distance,designation], where 'distance' is
%                  the distance to the supernova in kpc and 'designation'
%                  is the waveform type.  Known waveform types are 
%                  'A1B1G1','A1B1G2','A1B1G3','A1B1G4','A1B1G5','A1B2G1',
%                  'A1B2G2','A1B2G3','A1B2G4','A1B2G5','A1B3G1','A1B3G2',
%                  'A1B3G3','A1B3G4','A1B3G5','A2B1G1','A2B1G2','A2B1G3',
%                  'A2B1G4','A2B1G5','A2B2G1','A2B2G2','A2B2G3','A2B2G4',
%                  'A2B2G5','A2B3G1','A2B3G2','A2B3G3','A2B3G4','A2B3G5',
%                  'A2B4G2','A2B4G3','A2B4G4','A2B4G5','A2B5G4','A2B5G5',
%                  'A3B1G1','A3B1G2','A3B1G3','A3B1G4','A3B1G5','A3B2G1',
%                  'A3B2G2','A3B2G3','A3B2G4','A3B2G5','A3B3G1','A3B3G2',
%                  'A3B3G3','A3B3G4','A3B3G5','A3B4G2','A3B4G3','A3B4G4',
%                  'A3B4G5','A3B5G4','A3B5G5','A4B1G1','A4B1G2','A4B1G3',
%                  'A4B1G4','A4B1G5','A4B2G1','A4B2G2','A4B2G3','A4B2G4',
%                  'A4B2G5','A4B3G1','A4B3G2','A4B3G3','A4B3G4','A4B3G5',
%                  'A4B4G2','A4B4G3','A4B4G4','A4B4G5','A4B5G4','A4B5G5'.
%           The parameters may also be supplied as a single tilde-delimited
%           string.  MakeWaveform will convert this string to a cell or
%           double array using tildedelimstr2cellornum.
%   T       Scalar. Duration of the waveform.  Must have T*fs = integer.
%           Recommend T>=1 for astrophysical waveforms.
%   T0      Scalar.  Desired peak time of the waveform, as measured in hrss 
%           for (hp.^+hc.^2).^(0.5).  Must have 0<T0<T.  See "Notes" below 
%           for caveats.
%   fs      Scalar.  Sampling rate of the time series.  T*fs must be an 
%           integer.
%
%   t       Times at which the waveform is sampled, starting from zero.
%   hp      Waveform values in the plus polarization, in strain.  
%   hc      Waveform values in the cross polarization, in strain.
%
% The 'PropertyName','PropertyValue' pairs allow for optional arguments.
% Recognized pairs are:
%
%   'hrss',scalar 
%
%           Rescale waveform so that its root-sum-square (RSS) amplitude
%           has this value; i.e., 
%             \int dt hp(t)^2+hc(t)^2 = hrss^2
%           This over-rides the amplitude specified or implied by the 
%           waveform type-specific parameters.
%
%   'catalog', waveform_catalog
%
%           Search for waveform in the specified waveform catalog.  Useful
%           if you want to make many waveforms and reloading the catalog is
%           slow.
%
%   'catalogDirectory',string
%
%           Absolute path to the directory containing astrophysical
%           waveform catalogs (e.g., DFM, OB, or ZM).  If not supplied when
%           constructing such a waveform, xmakewaveform will attempt to
%           load the appropriate catalog from a hard-coded location.
%
% Notes:
%
% The waveform is interpolated so that the RSS-weighted center time,  
% [\int dt t (hp.^+hc.^2)].^(0.5), is at T0.  Note that for 2-polarization
% waveforms that include an inclination angle (eg, DS2P) different
% inclination angles will cause hp, hc to be weighted differently and
% therefore will produce different shifts.
%
% The cusp waveform (generated by P. Shawhan's cuspgw function) and the 
% astrophysical waveforms that are loaded from files are always 1 sec long
% with sampling rate 16384 Hz when created/loaded. These are resampled
% linearly and truncated or zero-padded to the desired length.  One should
% be cautious in using T<1 for these waveforms.
%
% $Id: xmakewaveform.m 1941 2007-11-14 05:00:57Z jrollins $

%----- Useful constants
pc = 3.08567756707701e+016;  %-- 1 parsec (m)
Mpc = 1e6*pc;  %-- 1 Mega-parsec (m)
c = 299792458;  %-- speed of light (m/s)

% ---- Assign default values to optional arguments.
hrss = -1;
Xcat = [];
%----- Hardwired waveform catalog storage directory.
filedir = ['/archive/home/mwas/lscsoft/src/matapps/src/searches/burst/coherent-network/waveforms' filesep];

% ---- Check for optional arguments.
if (nargin>5 && length(varargin))
    % ---- Make sure they are in pairs.
    if (length(varargin) == 2*floor(length(varargin)/2))
        % ---- Parse them.
        index = 1;
        while index<length(varargin)
            switch lower(varargin{index})
                case 'hrss'
                    hrss = varargin{index+1};
                case 'catalog'
                    Xcat = varargin{index+1};
                case 'catalogdirectory'
                    filedir = [varargin{index+1} filesep];
                otherwise
                    error(['Property value ' varargin{index} ' not recognized.']);
            end
            index = index + 2;
        end
    else
        error('Optional arguments must appear in pairs: property_name, property_value.');
    end
end

%----- Boolean: If true then waveform "type" specifies a waveform made 
%      outside of this function (loaded from file, or using another 
%      function) and which therefore may need to be re-sampled, shifted in 
%      time, zero-padded, or truncated.
pregen = 1;      %-- By default perform conditioning on everything.
pregen_fs = fs;  %-- sampling rate of pregenerated waveform.
pregen_T = T;    %-- duration of pregenerated waveform.

%----- Time in column vector.
t = [0:1/fs:T-1/fs]';

%----- If parameters are supplied as tilde-delimited string, then convert
%      to cell or double array.
if (ischar(parameters)) 
    parameters = tildedelimstr2numorcell(parameters); 
end

%----- Make desired waveform.
switch type
    
    case {'CG','cg'}

    %----- cosine-Gaussians
    h_peak = parameters(1);
    tau = parameters(2);
    f0 = parameters(3);
    hp = h_peak*cos(2*pi*(t-T0)*f0).*exp(-(t-T0).^2./tau.^2);
    hc = zeros(size(hp));

    %----- Turn off default interpolation (symmetric ad hoc waveform).
    pregen = 0;

    case 'cusp'

    %----- Cusp GWB made using function by Peter Shawhan.
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;
    hp = cuspgw(parameters(1));
    hp = hp(:);
    hc = zeros(size(hp));

    case 'chirplet'

    %----- Chirplet - Gaussian-modulated sinusoid with frequency changing
    %      linearly with time.  Put chirping cosine-Gaussian in plus 
    %      polarization, chirping sine-Gaussian in cross.
    h_rss = parameters(1);
    tau = parameters(2);
    f0 = parameters(3);
    alpha = parameters(4);
    delta = parameters(5);
    % chirpsign = parameters(6);
    % %----- Checks on parameters.
    % if (abs(chirpsign)~=1)
    %     error('Chirp sign must be +1 or -1.');
    % end
    % V = 4*pi*tau*bw;
    % if (V<1)
    %     error('Parameters unphysical.  Must have duration*bandwith>=0.5.');
    % end
    % alpha = chirpsign*(V^2-1)^0.5;
    h = h_rss*exp(...
            (-1+i*alpha)*(t-T0).^2./(4*tau.^2) ...
            +i*2*pi*(t-T0)*f0 ...
            +i*delta  ...
        )./(2*pi*tau^2).^(1/4);
    hp = real(h);
    hc = imag(h);

    %----- Turn off default interpolation (ad hoc waveform is designed to
    %      produce desired T0).
    pregen = 0;

    case {'DFM','dfm'}

    %----- Dimmelmeier-Font-Mueller supernova waveform.
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;

    %----- Load waveform structure, if not already loaded or supplied.
    if (isempty(Xcat))
        %disp('Loading Dimmelmeier-Font-Muller catalog.') 
        load([filedir 'DFMcatv6.mat'])
    end

    %----- Find specified waveform type.
    NAME = parameters{1,2};
    for k=1:length(DFMcat)
        if (strcmp(DFMcat(k).name,NAME))
            DFM_number = k;
        end
    end

    %----- Read waveform, which is defined at a range of 1kpc.
    %      Make sure h is same type of vector (column) as t.
    distance = 1e-3*parameters{1};
    h = 1./distance*DFMcat(DFM_number).hoft';
    hp = h(:);
    hc = zeros(size(hp));

    case {'DFMc','dfmc'}

    %----- Dimmelmeier-Font-Mueller supernova waveform.
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;

    %----- Load waveform structure, if not already loaded or supplied.
    if (isempty(Xcat))
        %disp('Loading Dimmelmeier-Font-Muller catalog.') 
        load([filedir '/DFMcatv6.mat'])
    end

    %----- Find specified waveform type.
    NAME = parameters{1,2};
    for k=1:length(DFMcat)
        if (strcmp(DFMcat(k).name,NAME))
            DFM_number = k;
        end
    end

    %----- Read waveform, which is defined at a range of 1kpc.
    %      Make sure h is same type of vector (column) as t.
    distance = 1e-3*parameters{1};
    h = 1./distance*DFMcat(DFM_number).hoft';
    % ---- Detrend.
    % ---- Get all nonzero samples of the waveform
    index = find(h ~= 0);
    index_start = index(1);
    index_end = index(end);
    % ---- In case there are zero values "inside"... 
    index = [index_start:index_end];
    hnz = h(index);
    % ---- Remove linear trend (set first and last points to zero)
    trend = hnz(end)*([0:length(hnz)-1]')/(length(hnz)-1); 
    ht = hnz - trend;
    trend = ht(1)*([length(hnz)-1:-1:0]')/(length(hnz)-1); 
    ht = ht - trend;
    h(index) = ht;
    %
    hp = h(:);
    hc = zeros(size(hp));

    case {'DS','ds'}

    %----- damped sinusoids
    h_peak = parameters(1);
    tau = parameters(2);
    f0 = parameters(3);
    hp = zeros(size(t));
    hc = zeros(size(t));
    k = find(t>=T0);
    hp(k) = h_peak*cos(2*pi*(t(k)-T0)*f0).*exp(-(t(k)-T0)/tau);
    hc(k) = h_peak*sin(2*pi*(t(k)-T0)*f0).*exp(-(t(k)-T0)/tau);

    %----- Force interpolation to get correct T0, since waveform is 
    %      asymmetric.
    pregen = 1;
    
    case {'DS2P','ds2p'}
    
    %----- Damped sinusoids, 2-polarization, with smooth turn-on.
    h_peak = parameters(1);
    tau = parameters(2);
    f0 = parameters(3);
    delta = parameters(4);
    ciota = parameters(5);

    %----- Smooth turn on: prepend time-reversed damped sinusoid with same
    %      parameters, except tau->tau/10.
    tau_vec = tau*ones(size(t));
    k = find(t<T0);
    tau_vec(k) = tau_vec(k)/10;
    hp = h_peak*0.5*(1+ciota^2)*cos(2*pi*(t-T0)*f0+delta).*exp(-abs(t-T0)./tau_vec);
    hc = h_peak*ciota*sin(2*pi*(t-T0)*f0+delta).*exp(-abs(t-T0)./tau_vec);

    %----- Force interpolation to get correct T0, since waveform is 
    %      asymmetric.
    pregen = 1;

    case {'G','g'}

    %----- Gaussians
    h_peak = parameters(1);
    tau = parameters(2);
    hp = h_peak*exp(-(t-T0).^2/tau.^2);
    hc = zeros(size(hp));
    
    %----- Turn off default interpolation (symmetric ad hoc waveform).
    pregen = 0;

%          'INSP'  10+10 solar-mass inspiral and merger waveform, at some 
%                  unknown distance (!).  See
%                  http://www.lsc-group.phys.uwm.edu/cgi-bin/bag-enote.pl?nb=burs3sim&page=3
%           INSP:  [].  Waveform has plus polarization only.  Must use T=1,
%                  fs=16384.  T0, Xcat ignored.
%     case {'INSP','insp'}
% 
%     %----- Kludged 10+10 M_o inspiral + merger waveform.  See
%     %      http://www.lsc-group.phys.uwm.edu/cgi-bin/bag-enote.pl?nb=burs3sim&page=3
%     %----- Turn off default interpolation.
%     pregen = 0;
%     % pregen_fs = 16384;
%     % pregen_T = 1;
% 
%     %----- Load waveform structure.  Distance is unknown.
%     filedir = 'C:\cygwin\home\psutton\lscsoft\matapps\src\searches\burst\coherent-network\waveforms\';
%     hp = load([filedir 'inspiral_10_10.txt']);
%     hp = hp(:);
%     hc = zeros(size(hp));
    
    case {'Lazarus','lazarus'}

    %----- Lazarus Project black-hole mergers

%     %----- Load waveform for 1 solar mass binary system.
%     %      column 1: t (s)
%     %      column 2: h_plus*dist (cm)
%     %      column 1: h_cross*dist (cm)
%     %      All three columns scale with system mass in solar masses.
%     m1 = load([filedir 'Lazarus.txt']);
%     %----- Rescale to get waveform for desired mass, distance.
%     mass = parameters(1);
%     dist = parameters(2);
%     ciota = parameters(3);
%     t_temp = m1(:,1)*mass;
%     hp_temp = m1(:,2)*mass/dist*(0.01/Mpc)*(0.5*(1+ciota^2));
%     hc_temp = m1(:,3)*mass/dist*(0.01/Mpc)*(ciota);
%     hp = hp_temp(:);
%     hc = hc_temp(:);
% 
%     %----- Need to inteprolate to desired sampling rate, etc.
%     pregen = 1;
%     pregen_T = length(t_temp)*(t_temp(end)-t_temp(1))/(length(t_temp)-1);
%     pregen_fs = length(t_temp)/pregen_T;

    % Use analytic approximation from Baker, Campanelli, Lousto, and 
    % Takahashi, PRD 65 124012 2002.

    % Input parameters
    mass_solarMasses = parameters(1);
    distance_megaParsecs = parameters(2);
    ciota = parameters(3);

    % Physical constants
    NewtonG = 6.67e-11;
    speedOfLight_MetersPerSecond = 299792458;
    megaParsec_Meters = 1e6*3.26*speedOfLight_MetersPerSecond*365.25*86400;
    solarMass_Kilograms = 1.99e30;
    solarMass_Seconds = NewtonG*solarMass_Kilograms/speedOfLight_MetersPerSecond^3;
    solarMass_Meters = NewtonG*solarMass_Kilograms/speedOfLight_MetersPerSecond^2;

    % Constants of analytic approximation.
    alpha0 = 0.0085;
    omega0 = 0.2894;
    tomega0 = 33.00;
    tomega1 = 67.05;
    sigma0 = 0.2192;
    tsigma0 = 32.54;
    tsigma1 = 62.92;
    omegaQNM = 0.55;
    sigmaQNM = -0.073;
    phi0 = -4.40;
    a0 = -6.31;

    % This range of t includes the entire range of envelope A > 1e-4*max(A).
    dt = 0.1;  % this is fine enough resolution (in units of system mass)
    tmin = 0;
    tmax = 200; % range over which analytic approximation has support.
    t = [tmin:dt:tmax]';

    % Construct phase phi(t) piecewise in time.
    phi = zeros(size(t));
    k = find(t<tomega0);
    phi(k) = phi0 + omega0*(t(k)-tomega0) + alpha0/2*(t(k)-tomega0).^2;
    k = find( (t>=tomega0) & (t<tomega1) );
    phi(k) = phi0 + omega0*(t(k)-tomega0) + (omega0-omegaQNM)/(tomega0-tomega1)/2*(t(k)-tomega0).^2;
    k = find(t>=tomega1);
    phi(k) = phi0 + 1/2*(omega0+omegaQNM)*(tomega1-tomega0) + omegaQNM*(t(k)-tomega1);

    % Construct amplitude A(t) piecewise in time.
    lnA = zeros(size(t));
    k = find(t<tsigma0);
    lnA(k) = a0 + sigma0*(t(k)-tsigma0);
    k = find( (t>=tsigma0) & (t<tsigma1) );
    lnA(k) = a0 + sigma0*(t(k)-tsigma0) + 1/2*(sigma0-sigmaQNM)/(tsigma0-tsigma1).*(t(k)-tsigma0).^2;
    k = find(t>=tsigma1);
    lnA(k) = a0 + 1/2*(sigma0+sigmaQNM)*(tsigma1-tsigma0) + sigmaQNM*(t(k)-tsigma1);
    A = exp(lnA);

    % Construct r*Psi_4 variable:
    rPsi4 = A.*exp(-i*phi);

    % Integrate twice in time to get hplus, hcross.
    ddhp = real(2*rPsi4); 
    ddhc = imag(2*rPsi4); 
    dhp = dt*real(cumsum(2*rPsi4)); 
    dhc = dt*imag(cumsum(2*rPsi4)); 
    hp = dt^2*real(cumsum(cumsum(2*rPsi4)));
    hc = dt^2*imag(cumsum(cumsum(2*rPsi4)));
    % remove linear trend in hp:
    tint = t(end)-hp(end)/dhp(end);
    m = dhp(end);
    b = hp(end) - dhp(end)*t(end);
    hp_trend = zeros(size(t));
    k = find(t>tint);
    hp_trend(k) = m*t(k)+b;
    hp = hp - hp_trend;
    % remove linear trend in hc:
    tint = t(end)-hc(end)/dhc(end);
    m = dhc(end);
    b = hc(end) - dhc(end)*t(end);
    hc_trend = zeros(size(t));
    k = find(t>tint);
    hc_trend(k) = m*t(k)+b;
    hc = hc - hc_trend;

    % Rescale to seconds and megaParsecs.
    t = t * mass_solarMasses * solarMass_Seconds;
    hp = hp * ( mass_solarMasses * solarMass_Meters ) / ( distance_megaParsecs * megaParsec_Meters );
    hc = hc * ( mass_solarMasses * solarMass_Meters ) / ( distance_megaParsecs * megaParsec_Meters );

    % rescale waveforms according to inclination angle.
    hp = hp*(0.5*(1+ciota^2));
    hc = hc*ciota;

    %----- Need to inteprolate to desired sampling rate, etc.
    pregen = 1;
    pregen_T = length(t) * dt * mass_solarMasses * solarMass_Seconds;
    pregen_fs = length(t)/pregen_T;

    case {'OB','ob'}

    %----- Ott-Burrows-et_al supernova waveform (pregenerated).
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;

    %----- Load waveform structure, if not already loaded or supplied.
    if (isempty(Xcat))
        %disp('Loading Ott-Burrows catalog.') 
        load([filedir 'OBcatv6.mat'])
    end

    %----- Find specified waveform type.
    NAME = parameters{1,2};
    for k=1:length(OBcat)
        if (strcmp(OBcat(k).name,NAME))
            OB_number = k;
        end
    end

    %----- Read waveform, which is defined at a range of 10kpc.
    %      Make sure h is same type of vector (column) as t.
    distance = 0.1*parameters{1};
    h = 1./distance*OBcat(OB_number).hoft';
    hp = h(:);
    hc = zeros(size(hp));

    case 'onecyclesine'

    %----- one cycle of a sine wave.
    h_peak = parameters(1);
    f0 = parameters(2);
    hp = h_peak*sin(2*pi*(t-T0)*f0);
    k = find(abs(t-T0)>=1/abs(2*f0));
    hp(k) = 0;
    hc = zeros(size(hp));
    
    %----- Turn off default interpolation (symmetric ad hoc waveform).
    pregen = 0;

    case {'SG','sg'}

    %----- sine-Gaussians
    h_peak = parameters(1);
    tau = parameters(2);
    f0 = parameters(3);
    hp = h_peak*sin(2*pi*(t-T0)*f0).*exp(-(t-T0).^2./tau.^2);
    hc = zeros(size(hp));
    
    %----- Turn off default interpolation (symmetric ad hoc waveform).
    pregen = 0;

    case 'whistle'
    
    duration = parameters(1);
    fstart = parameters(2);
    fend = parameters(3);
    hp = whistlegw(duration,fstart,fend);
    hp = hp(:);
    hc = zeros(size(hp));

    %----- Turn on interpolation (waveform generated by outside code).
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;

    case {'WNB','wnb'}

    %----- Turn off default interpolation.
    %      (Might actually be useful here, but don't want to alter
    %      frequency content of noise by low-pass filtering).
    pregen = 0;

    %----- Gaussian-modulated noise burst, white over specified band.
    h_rss = parameters(1);
    fc = parameters(2);
    df = parameters(3);
    dt = parameters(4);
    %----- Gaussian envelope
    env = exp(-(t-T0).^2/2/dt.^2);
    %----- Band-limited noise (independent for each polarization)
    x = BLWNB(max(fc-df,0),2*df,T,fs);
    x = x';
    hp = env.*x;
    hp = hp*h_rss/(hp'*hp/fs).^0.5;
    x = BLWNB(max(fc-df,0),2*df,T,fs);
    x = x';
    hc = env.*x;
    hc = hc*h_rss/(hc'*hc/fs).^0.5;    

    case 'zero'

    %----- Trivial (null) signal.  Useful to avoid scripts crashing.
    hp = zeros(size(t));
    hc = zeros(size(t));
    
    %----- Turn off default interpolation.
    pregen = 0;

    case {'ZM','zm'}

    %----- Zwerger-Mueller supernova waveform (pregenerated).
    pregen = 1;
    pregen_fs = 16384;
    pregen_T = 1;

    %----- Load waveform structure, if not already loaded or supplied.
    if (isempty(Xcat))
        %disp('Loading Zwerger-Muller catalog.') 
        load([filedir 'ZMcatv6.mat'])
    end

    %----- Find specified waveform type.
    NAME = parameters{1,2};
    for k=1:length(ZMcat)
        if (strcmp(ZMcat(k).name,NAME))
            ZM_number = k;
        end
    end

    %----- Read waveform, which is defined at a range of 1Mpc.
    %      Make sure h is same type of vector (column) as t.
    distance = 1e-3*parameters{1};
    h = 1./distance*ZMcat(ZM_number).hoft';
    hp = h(:);
    hc = zeros(size(hp));
    
    otherwise
    
    error(['Waveform type ' type ' not recognized.'])
    
end

%----- Resample, time-shift, truncate, and/or zero-pad pregenerated
%      waveform as needed.
if (pregen)

    %----- Time vector for pregenerated waveform.
    pregen_t = [0:1/pregen_fs:pregen_T-1/pregen_fs]';
    %----- Measure peak/characteristic time of pregenerated waveform.
    [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = xoptimalsnr( ... 
        (hp.^2+hc.^2).^0.5,pregen_t(1),pregen_fs,[],[],[], ...
        1/pregen_T,0.5*pregen_fs-1/pregen_T ...
    );
    %----- Shift pregenerated time vector by required amount.
    pregen_t = pregen_t + (T0 - Tchar);
    %----- Interpolate to correct sampling rate and peak time.
    %----- Time in column vector.
    t = [0:1/fs:T-1/fs]';
    %----- Find desired times which overlap pregen_t.
    k = find(t>=pregen_t(1) & t<=pregen_t(end));
    %----- Interpolate, with zero padding if pregenerated waveform is too
    %      short.  (Truncation of long waveforms handled automatically by
    %      specifying vector "t".) 
    hp_interp = zeros(size(t));
    hc_interp = zeros(size(t));
    hp_interp(k) = interp1(pregen_t,hp,t(k));
    hc_interp(k) = interp1(pregen_t,hc,t(k));
    hp = hp_interp;
    hc = hc_interp;
    
end

%----- Reset hrss amplitude to specified value, if any. 
%      For 2-polarization waveforms use hp^2+hc^2.
%      Note: If time-shifting moves pregen waveform out of interval of
%      interest then this will screw up the amplitude!
if ( (nargin>=6) && (isempty(hrss)==0) && (hrss>0) )
    norm = ((hp'*hp+hc'*hc)/fs)^0.5;
    hp = hp*hrss/norm;
    hc = hc*hrss/norm;
end

%----- Done.
return
