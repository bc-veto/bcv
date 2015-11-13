function moments = wmoment(triggers, momField, momNum, weightField, ...
                           normalization)
% WMOMENT Calculate a particular (central) moment from a set of triggers
%
% WMOMENT Evaluates an arbitrary (central) moment integral for a particular field,
% weighted by a given weighting field. The normalization may be specified.
%
% usage:  moments = wmoment(triggers, momField, momNum, weightField, ...
%                           normalization);
%
%   triggers          cell array of significant tile properties.
%   momField          string name of field specifying which moment to evaluate.
%   momNum            scalar specifying the power n in <(x-mu)^n> for the integral
%   weightField       (optional) scalar name for the weighting function 
%                     for the integral. If not specified, defaults to
%                     normalizeEnergy.
%   normalization     (optional) vector for the normalization of the
%                     integrals. If not specified, defaults to the integral of
%                     the weighting function alone, for each channel.
%
% Note that since the transform basis tiles each have dt*df = 1, all of the
% integrations here omit the t-f area element for convenience.
%
% For n==1, the moment is not a central moment.
%
% See also WTRANSFORM, WTHRESHOLD, and WSELECT.

% Leo C. Stein
% lstein@ligo.caltech.edu

% $Id: wmoment.m 2263 2009-08-24 19:54:40Z jrollins $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 5, nargin));

% apply default arguments
if (nargin < 3) || isempty(momNum),
  momNum = 1;
end
if (nargin < 4) || isempty(weightField),
  weightField = 'normalizedEnergy';
end
if (nargin < 5) || isempty(normalization),
  normalization = []; % taken care of below.
end


% force cell arrays
triggers = wmat2cell(triggers);

% force one dimensional cell arrays
triggers = triggers(:);

% determine number of channels
numberOfChannels = length(triggers);

% preallocate space for moments
moments = zeros(numberOfChannels, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate triggers structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(triggers{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
  if ~isfield(triggers{channelNumber},momField),
    error(['no field named "' momField '" in channel ' num2str(channelNumber)]);
  end
  if ~isfield(triggers{channelNumber},weightField),
    error(['no field named "' momField '" in channel ' num2str(channelNumber)]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            switch for first moment                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if momNum ~= 1,
  mus = wmoment(triggers, momField, 1, weightField, normalization);
else
  mus = zeros(numberOfChannels, 1);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  curChan = triggers{channelNumber};
  mu = mus(channelNumber);
  
  if isempty(normalization),
    curNorm = sum(curChan.(weightField));
  else
    curNorm = normalization(channelNumber);
  end
  
  moments(channelNumber) = sum((curChan.(momField) - mu) .^ momNum .* ...
                               curChan.(weightField)) / curNorm;
  
% end loop over channels
end