function event = ...
        wcoherentclusteranalysis(data, tiling, coefficients, event, ...
                                 triggers, candidateChannel, ...
                                 blockStartTime, energyTypes, debugLevel)
% This function essentially takes output from the bayesian code and uses it 
% to calculate the energies for the pixels in the selected cluster.
%
% Also acts as a stable interface buffer for wcoherentanalysis.
% It packages all needed data into neat structures before calling the main
% function.
% At the moment, it also handles the loop over the cluster pixels.
% This will either be incorporated into wcoherentanalysis, or handled in a
% better way - i.e. check pixel.duration and pull out the data associated with
% those pixels in a single tfmap, rather than creating anew each time.

% Mark Edwards <Mark.Edwards@astro.cf.ac.uk>

% $Id:$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Process input arguments                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(triggers);

% get channel names, which are available in the full trigger structure
for channelNumber = 1 : numberOfChannels,
  channelNames{channelNumber} = triggers{channelNumber}.channelName;
end

% now that we have the channel names, we only need the triggers from the
% candidate channel
triggers = triggers{candidateChannel};

% Add the energy fields to the event structure - allows the generation of a
% running sum in the pixel cluster loop
for energyNumber = 1:length(energyTypes)
  event.([energyTypes{energyNumber} 'Energy']) = 0;
  event.([energyTypes{energyNumber} 'IncoherentEnergy']) = 0;
end

% Extract sky position for event - calculated by the Bayesian code (wrap around 
% if 0>theta>pi or 0>phi>2*pi)
skyPosition = [mod(event.modeTheta,pi), mod(event.modePhi,2*pi)];

% Get indices of the triggers that belong to the cluster
pixelIndex = find(event.clusterNumber == triggers.clusterNumber);

% Get number of events in this cluster - at the moment, event does not seem to
% return the corrent cluster size, so do not rely on it - use pixelIndex too
noPixels = min(event.size, length(pixelIndex));

% Loop around pixels in cluster, adding the calculated energies to the energy 
% field in the event structure 

for eachPixel = 1:noPixels
  % Need to retrieve all the information needed to run wcoherentanalysis.
  % Such as frequency, duration etc. Get this from triggers struct

  % Also, pack data for a call to the function - use my own structure to provide
  % a stable interface to the function
  pixels(eachPixel).q         = triggers.q(pixelIndex(eachPixel));
  pixels(eachPixel).frequency = triggers.frequency(pixelIndex(eachPixel));
  pixels(eachPixel).time      = triggers.time(pixelIndex(eachPixel));
  pixels(eachPixel).duration  = triggers.duration(pixelIndex(eachPixel));
  
  properties.blockStartTime      = blockStartTime;
  properties.theta               = skyPosition(1);
  properties.phi                 = skyPosition(2);
  properties.sampleFrequency     = tiling.sampleFrequency;
  properties.blockDuration       = tiling.duration;
  properties.candidateChannel    = candidateChannel;
  properties.channelNames        = channelNames;
  
  % Calculate noise in our frequency band
  % Get indices for the plane and row corresponding to the pixel
  best_q = find(tiling.qs == triggers.q(pixelIndex(eachPixel)));
  best_row = find (tiling.planes{best_q}.frequencies == triggers.frequency(pixelIndex(eachPixel)));
  row = tiling.planes{best_q}.rows{best_row};

  wSw = zeros(numberOfChannels,1);
  
  % Generate PSD, or in this case 1/PSD.
  for channelNumber = 1:numberOfChannels
    % Normalize window
    row.window = row.window / sqrt(sum(row.window.^2));
    
     % Calculate wSw for each detector in our frequency band 
     wSw(channelNumber) = sum((row.window .* coefficients{channelNumber}(row.dataIndices)).^2)/2;
  end

  % Call coherent analysis with required data - returns the event structure with 
  % the energy fields populated with the cluster energies, and a pixelEnergy 
  % structure which has the calculated energies of each of the pixels
  [event, pixelEnergy(eachPixel)] = wcoherentanalysis(wSw, event, data, ...
            coefficients, energyTypes, pixels(eachPixel), properties);

end
