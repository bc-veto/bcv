function wwriteskymap(skymap, filename)
% WWRITESKYMAP Write a skymap to file
%
% wwriteskymap(skymap, filename)
%
% input:
%   skymap     a skymap
%   filename   name of file to be written
%
% format:
%   phi theta logProbability

if (any(any(isnan(skymap))))
  error('NaNs in the skymap!');
end

% return only the top 100 square degrees
n = ceil(length(skymap)*100/(360*180*2/pi));

% convert from log probability to probability, normalized to unity on sky

% exponentiate, scaling so largest value is one (to avoid overflow)
a = exp(skymap(:, 3) - skymap(1, 3));
% normalize so skymap sums to unity
skymap(:, 3) = a / sum(a);

% write the file out
fid = fopen(filename, 'w');
for k=1:n,
  fprintf(fid, '%.10e %.10e %.10e\n', skymap(k,1), skymap(k,2), skymap(k,3));
end
fclose(fid);
