function [logSkymap, logGlitches] = wposteriorsgs(channelNames, rate, wSw, xSw, directions, taus)
%
% Performs Bayesian followup for various signal/glitch sizes

% first, we normalize everything so that max(wSw) == 1

a = max(wSw);

wSw = wSw / a;

for i = 1:length(xSw)
    xSw{i} = xSw{i} / sqrt(a);
end

% now we are set up to analyze for a signal at SNR 1

sigmaCount = 6;
sigmaFactor = 10;

% 1.00 3.16 10.0 31.6 100 316

for j = 1:sigmaCount

    % generate the posteriors

    if (nargout == 2)
    
        [b, c] = wposteriorsg(channelNames, rate, wSw, xSw, directions, taus);
    
    else
        
        b = wposteriorsg(channelNames, rate, wSw, xSw, directions, taus);
        c = 0;    
    end
    
    %figure(j + 1)
    %imagesc(reshape(b, 361, 181)')
    %figure(1)

    % accumulate the posteriors
    
    
    if (j == 1)
        
        % if this is the first time, initialize the outputs
        
        logSkymap = b;
        logGlitches = c;
               
    else
        
        % accumulate with the existing values
        
        logSkymap = logsumexp(logSkymap, b);
        logGlitches = logsumexp(logGlitches, c);
        
    end
    
    % scale up w for the next pass
    
    wSw = wSw * sigmaFactor;
    for i = 1:length(xSw)
        xSw{i} = xSw{i} * sqrt(sigmaFactor);
    end
    
end

% normalize the results

logSkymap = logSkymap - log(sigmaCount);
logGlitches = logGlitches - log(sigmaCount);
