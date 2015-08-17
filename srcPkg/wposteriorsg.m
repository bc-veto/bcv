function [logSkymap, logGlitches] = wposteriorsg(channelNames, rate, wSw, xSw, directions, taus)
%
% Performs Bayesian followup using wsposteriors to compute skymap and
% glitch probabilities
%
% first, the coherent skymap

logSkymap = wposteriors(channelNames, rate, wSw, xSw, directions, taus);

% next, the incoherent skymap for each glitch

if nargout == 2

for i = 1:length(channelNames)
    
    gtaus = [ taus(i * 2 - 1), taus(i * 2), taus(end - 2), taus(end - 1), taus(end) ];
    
    a = wposteriors(channelNames(i), rate, wSw(i), xSw(i), directions, gtaus);
    
    % only store the glitch probability
    
    logGlitches(i) = log(mean(exp(a - max(a)))) + max(a);
    
end

end

