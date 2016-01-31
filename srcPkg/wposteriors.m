% logSkymap = wposteriors(channelNames, rate, wSw, xSw, directions, taus)
%
% Performs bayesian followup for a timeseries
%
% taus = (rate/4):(rate*3/4);
% channelNames = { 'H1:foo', 'L1:bar', 'V1:snafu' };
% rate = 4096;
% wSw = [ 1, 1, 1 ];
% xSw = { randn(1, rate), randn(1, rate), randn(1, rate) };
% [theta, phi] = meshgrid(0:0.1:pi, -pi:0.1:pi);
% directions = [ theta(:)' ; phi(:)' ];
% taus = (rate/4):(rate*3/4);




% [skymap, posteriors] = wposteriors([duration sampleRate analysisRate minFrequency maxFrequency minTime maxTime], psd, x)
%
% Performs Bayesian followup.  Under devlopment.  
%
% Test invocation:
%
% duration = 64;
% sampleRate = 4096;
% analysisRate = 4096;
% minFrequency = 100;
% maxFrequency = 300;
% minTime = 31.9;
% maxTime = 32.1;
% parameters = [duration sampleRate analysisRate minFrequency maxFrequency minTime maxTime];
% m = sampleRate * duration / 2;
% psd  = {ones([m 1]), ones([m 1]), ones([m 1])};
% data = {randn([m 1]) + j * randn([m 1]), ...
%         randn([m 1]) + j * randn([m 1]), ...
%         randn([m 1]) + j * randn([m 1]) };
% [skymap, posteriors] = wposteriors(parameters, psd, data);
% lognoise  = 1
% logsignal = posteriors(1)
% logglitch = posteriors(2)
% modetheta = posteriors(3)
% modephi   = posteriors(4)
% properties = reshape(posteriors(3:end), 8, []);
% bands.minFrequency = properties(1,:);
% bands.maxFrequency = properties(2,:);
% bands.site{1}.minTime = properties(3, :);
% bands.site{1}.maxTime = properties(4, :);
% bands.site{2}.minTime = properties(5, :);
% bands.site{2}.maxTime = properties(6, :);
% bands.site{3}.minTime = properties(7, :);
% bands.site{3}.maxTime = properties(8, :);
% imagesc(skymap)
%
% where 
%    duration, minTime, maxTime are in seconds
%    sampleRate, analysisRate, minFrequency, maxFrequency are in Hz
%    psd and data are one-sided

error('wposteriors.m was executed; please compile wposteriors.cpp into a mex-file');

