function [injectionData, gps_s, gps_ns, phi, theta, psi] = winject( ...
    startTime,blockTime,channel,sampleRate,injectionFileName, ...
    rescaleByAntennaResponse,injectionNumber,varargin)
% winject: Create simulated data streams containing simulated signals 
% for the specified network of gravitational-wave detectors.
%
% Copied from XINJECTSIGNAL
%
%  [injectionData, gps_s, gps_ns, phi, theta, psi] = xinject( ...
%    startTime,blockTime,channel,sampleRate,injectionFileName, ...
%    rescaleByAntennaResponse,injectionNumber)
%
%  startTime    Scalar.  Start time of the data being analysed.
%  blockTime    Scalar.  Duration (s) of the data being analysed.
%  channel      Cell array of strings.  channel{j} holds the name of the
%               jth detector channel being analysed.  The detector name is
%               inferred from the first character of the channel name, and
%               must be one of the names recognized by LoadDetectorData.
%  sampleRate   length(channel)x1 vector.  Sample rates corresponding to
%               channels.
%  injectionFileName
%               String.  Name of file specifying simulated signals to be
%               injected.
%  rescaleByAntennaResponse
%               Scalar.  Optional.  If nonzero then rescale each injection 
%               amplitude by 1/(sum over detectors of Fp^2)^0.5, where Fp
%               is calculated for each detector using the sky angles of the
%               injection into the first detector.  This rescaling is
%               useful for injecting linearly polarized GWBs into identical
%               detectors so that the sum-squared SNR is held fixed.
%               WARNING: It is probably not a good idea to use
%               this rescaling for: (1) glitches simulated using different
%               sky positions for each detector; (2) waveforms with 2
%               polarizations; (3) detectors that are not identical.
%  injectionNumber
%               Array (positive integer).  Optional.  If given then only
%               the injections specified on rows "injectionNumber" of the
%               injection file are injected.  This over-rides the default
%               behaviour, which is to inject all injections in the
%               interval [statTime,startTime+blockTime].
%
%  doNotInject  Logical. Optional. If set to true, waveform is not injected
%               only parameters are calculated
% 
%  injectionData
%               length(channel)x1 cell array of column vectors.
%               injectionData{j} is a column vector holding the signal data
%               for the jth detector.
%  gps_s        Vector of peak time (GPS seconds, integer part) of each 
%               injection at the center of the Earth.   
%  gps_ns       Vector of peak time (GPS seconds, nanosecond part) of each
%               injection at the center of the Earth.   
%  phi          Vector of azimuthal sky coordinate of each injection
%               (radians, Earth-fixed coordinates).   
%  theta        Vector of polar sky coordinate of each injection 
%               (radians Earth-fixed coordinates).   
%  psi          Vector of polarization angle of each injection  
%               (radians Earth-fixed coordinates).   
%
% The coordinate system is explained in ComputeAntennaResponse. It is
% assumed that the sky positions in the injection file are supplied in this
% coordinate system.  For glitch injections the reported times and angles
% are those for the injection into the first detector.
%
% WARNING: This simulation engine is designed for signals of duration less
% than approximately one second.  For longer signals the code will have to
% be modified slightly.
%
%  initial version: Patrick J. Sutton 2005.07.10
%
% $Id: xinjectsignal.m 1941 2007-11-14 05:00:57Z jrollins $


% ---- Check for optional arguments.
doNotInject=0;
catalogDirectory='';
if (nargin>7 && length(varargin))
    % ---- Make sure they are in pairs.
    if (length(varargin) == 2*floor(length(varargin)/2))
        % ---- Parse them.
        index = 1;
        while index<length(varargin)
            switch lower(varargin{index})
                case 'donotinject'
                    doNotInject = varargin{index+1};
                case 'catalogdirectory'
                    catalogDirectory = varargin{index+1};
                otherwise
                    error(['Property value ' varargin{index} ' not recognized.']);
            end
            index = index + 2;
        end
    else
        error('Optional arguments must appear in pairs: property_name, property_value.');
    end
end
if (nargin<7)
     injectionNumber = 0;
end
if (nargin<6)
     rescaleByAntennaResponse = 0;
end

%----- Speed of light (m/s).
c = 299792458;

%----- Parse channel list and load info on corresponding detectors.
nDetector = length(channel);
for j=1:nDetector
    detector{j,1} =  LoadDetectorData(channel{j}(1));
end

%----- Read injection file and extract parameters of injections during
%      time interval of interest.
injectionParameters = readinjectionfile(startTime,blockTime,injectionFileName);

%----- If specific injections are specified then keep only those injections.
if (injectionNumber)
    temp_injectionParameters = injectionParameters(injectionNumber);
    clear injectionParameters
    injectionParameters = temp_injectionParameters;
end

%----- Number of injections.
nInjection = size(injectionParameters,1);

%----- Separate each line (a string) into separate elements in a cell 
%      array.  For actual GWB injections (instead of glitches) elements for
%      jParam > 7 may be empty.
for jParam=1:(7*nDetector)
    for iCell=1:nInjection
        [currentParameters{jParam}{iCell} injectionParameters{iCell}] = ...
            strtok(injectionParameters{iCell});
    end
end

%----- Compute sum over detectors of Fp^2 for each injection, if desired.
if (rescaleByAntennaResponse)
    totalInjectedPower = zeros(nInjection,1);
    %----- Loop over simulated signals and record peak time, sky angles.
    for jInjection=1:nInjection
        %----- Use angles for first detector.
        phi    = str2num(currentParameters{3}{jInjection});
        theta  = str2num(currentParameters{4}{jInjection});
        psi    = str2num(currentParameters{5}{jInjection});
        %----- Compute Fp for each detector.
        for jDetector=1:nDetector;
            Fp = ComputeAntennaResponse(phi,theta,psi,detector{jDetector}.d);
            totalInjectedPower(jInjection) = ...
                totalInjectedPower(jInjection) + Fp^2;
        end
    end
end

%----- Loop over detectors and construct simulated signals. 
for jDetector=1:nDetector;
    %----- Make timeseries data for given detector.
    injectionData{jDetector,1} = zeros(blockTime*sampleRate(jDetector),1);
    %----- Loop over simulated signals.
    for jInjection=1:nInjection
        %---------- Extract parameters for current injection.
        if (isempty(currentParameters{(jDetector-1)*7+1}{jInjection}))
            %----- This injection is a GWB; use same params for all
            %      detectors.
            gps_s  = str2num(currentParameters{1}{jInjection});
            gps_ns = str2num(currentParameters{2}{jInjection});
            phi    = str2num(currentParameters{3}{jInjection});
            theta  = str2num(currentParameters{4}{jInjection});
            psi    = str2num(currentParameters{5}{jInjection});
            GWB_type = currentParameters{6}{jInjection};
            GWB_params = currentParameters{7}{jInjection};
        else
            %----- There are separate injection parameters for this
            %      detector.  Use them.
            gps_s  = str2num(currentParameters{(jDetector-1)*7+1}{jInjection});
            gps_ns = str2num(currentParameters{(jDetector-1)*7+2}{jInjection});
            phi    = str2num(currentParameters{(jDetector-1)*7+3}{jInjection});
            theta  = str2num(currentParameters{(jDetector-1)*7+4}{jInjection});
            psi    = str2num(currentParameters{(jDetector-1)*7+5}{jInjection});
            GWB_type = currentParameters{(jDetector-1)*7+6}{jInjection};
            GWB_params = currentParameters{(jDetector-1)*7+7}{jInjection};
        end

        %----- Unit vector pointing towards source.
        omega = [sin(theta).*cos(phi) ; sin(theta).*sin(phi) ; cos(theta)];

        %----- Time delay for incoming signal wrt center of Earth 
        %      (t_det - t_coe).
        delay = -detector{jDetector}.V'*omega/c;

        %----- Compute antenna response, including polarization angle.
        [Fp Fc] = ComputeAntennaResponse(phi,theta,psi,detector{jDetector}.d);

        %----- Rescale amplitude by sum_over_detectors of Fp^2, if desired:
        if (rescaleByAntennaResponse)
            Fp = Fp / sqrt(totalInjectedPower(jInjection));
            Fc = Fc / sqrt(totalInjectedPower(jInjection));
        end

        %----- KLUDGE: Make injection data in 3-second snippet.  This part
        %      will have to be changed for signals with durations greater
        %      than about 1 second.
        snippetDuration = 3;
        snippetPad = 1;
        %----- Peak Time of signal relative to startTime.
        peakTime = gps_s+1e-9*gps_ns+delay;  
        %----- Waveform plus and cross polarizations, in wave frame.
        if(~doNotInject)
            if(~isempty(catalogDirectory))
                [t,hp,hc] = xmakewaveform(GWB_type,GWB_params,snippetDuration, ...
                    snippetPad+peakTime-floor(peakTime),sampleRate(jDetector),...
                    'catalogDirectory',catalogDirectory);
            else
                [t,hp,hc] = xmakewaveform(GWB_type,GWB_params,snippetDuration, ...
                    snippetPad+peakTime-floor(peakTime),sampleRate(jDetector));
            end
            %----- Reset t to GPS time.
            t = t + floor(peakTime)-snippetPad;
            %----- Sample indices.
            injectionSamples = [1:length(t)]';
            %----- Drop samples which fall outside the desired time interval.
            k = find( (t<startTime) | (t>=startTime+blockTime) );
            if (~isempty(k))
                injectionSamples(k) = [];
                t(k) = [];
                hp(k) = [];
                hc(k) = [];
            end
            %----- Make sure injectionSamples is not empty so we don't crash
            if(isempty(injectionSamples))
                continue;
            end
            %----- Combine with antenna response to give signal.
            injectionSamples = injectionSamples ...
                + round((floor(peakTime)-snippetPad-startTime)*sampleRate(jDetector));
            injectionData{jDetector,1}(injectionSamples) = ...
                injectionData{jDetector,1}(injectionSamples) ...
                + Fp*hp + Fc*hc;
        end
    end
end

%----- Record for output the peak time and sky angles for each of the
%      injections.  For glitches record only the parameters for the first
%      detector.
%----- Prepare storage.  Re-use variable names....
gps_s = zeros(nInjection,1);
gps_ns = zeros(nInjection,1);
phi = zeros(nInjection,1);
theta = zeros(nInjection,1);
psi = zeros(nInjection,1);
%----- Loop over simulated signals and record peak time, sky angles.
for jInjection=1:nInjection
    gps_s(jInjection)  = str2num(currentParameters{1}{jInjection});
    gps_ns(jInjection) = str2num(currentParameters{2}{jInjection});
    phi(jInjection)    = str2num(currentParameters{3}{jInjection});
    theta(jInjection)  = str2num(currentParameters{4}{jInjection});
    psi(jInjection)    = str2num(currentParameters{5}{jInjection});
end

%----- Done
return
