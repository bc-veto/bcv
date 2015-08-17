function [theta, phi] = wsinusoidalprojection(n)
% [THETA, PHI] = WSINUSOIDALPROJECTION(N)
%
% Returns a set of coordinates theta(i) phi(i) describing a sinusoidal
% projection with N unique values of theta

% values of theta evenly spaced
a = pi * ((1:n) - .5) / n;

theta = [];
phi = [];

% for each value of theta (line of latitude)
for u = a
    
    % the length of the line of latitude, rounded to an integer number of
    % pixels
    b = round(2 * n * sin(u));
    
    % the phi coordinates of that many pixels
    c = 2 * pi * ((1:b) - .5) / b;
    
    % pack the coordinates of the new pixels into the theta and phi vectors
    % (theta value is the same for all on this line of latitude)
    theta = [ theta, repmat(u, size(c)) ];
    phi = [ phi, c ];
    
end

theta = theta(:);
phi = phi(:);

