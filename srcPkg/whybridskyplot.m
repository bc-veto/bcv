function whybridskyplot(skymap)
% WHYBRIDGSKYPLOT(SKYMAP) Plot a bayesian sky map
%
% Plots a skymap of the format [ THETA, PHI, LOGP ] by rasterizing into
% a bitmap whose size is guessed by the number or rows (pixels)
%
% Authors:
% Antony Searle <antony.searle@anu.edu.au>

% $Id: $

% sort the plotted quantity so that the higher values are rasterized later
% (overwriting lower values that have been written to the same pixel)

[logp, i] = sort(skymap(:,3));

% extract (sorted) theta and phi

theta = skymap(i, 1);
phi = skymap(i, 2);

% estimate the size of the bitmap to use based on the number of samples
max_v = ceil(sqrt(length(theta)/2));
max_u = max_v * 2;
%max_v = 360;
%max_u = 720;


% work out the (u, v) coordinates of each sample
u = 1 + floor(phi * max_u / (2 * pi));
v = 1 + floor(theta * max_v / pi);

% guard against values outside the bitmap
u(u < 1) = 1;
u(u > max_u) = max_u;
v(v < 1) = 1;
v(v > max_v) = max_v;

% make a zero bitmap flattened to one dimension
c = zeros(1, max_u * max_v) + min(logp);

% rasterize the plotted quantity into the bitmap
c(1 + (v - 1) + max_v * (u - 1)) = logp;

% unflatten the bitmap into two dimensions
c = reshape(c, max_v, max_u);

% compute the pixel center coordinates to tell IMAGESC how to scale the
% bitmap
x = 2 * pi * ((1:max_u) - .5) / max_u;
y = pi * ((1:max_v) - .5) / max_v;

% display the bitmap
imagesc(x, y, c);

% show the scale
colorbar;

% label the axes
xlabel('phi');
ylabel('theta');
