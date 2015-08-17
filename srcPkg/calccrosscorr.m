function [r, rMax, xCorr] = calccrosscorr(u, v)
%
% CALCCROSSCORR - calculate the linear correlation coefficient. 
%
% Usage: [r, rMax, xCorr] = calccrosscorr(u, v)
%
% u, v  : vectors - real or complex. u and v must be of the same size
% r     : linear correlation coefficient 
% rMax  : r maximised over time-shifs between the vectors
% xCorr : correlation function between u and v
%
%  r = Real[<u, v*> / (sqrt(<u, u*>) * sqrt(<v, v*>))]
%
% P. Ajith, 21.09.06
% 
% $Id: calccrosscorr.m 275 2009-10-30 05:11:22Z ajith $

% Check that the size of vectors u and v are the same.
if size(u) ~= size(v)
    error('Size of u and v must be the same\n');
end

% Calculate the linear cross correlation coefficient
uDotU  = sum(u.*conj(u));
vDotV  = sum(v.*conj(v));
uConjV = u.*conj(v)./(sqrt(uDotU).*sqrt(vDotV));
r      = real(sum(uConjV));

% calculate the correlation function
xCorr = ifft(uConjV)*length(u);
rMax = max(real(xCorr));
rMin = min(real(xCorr));

if abs(rMin) > abs(rMax)
    rMax = rMin;
end
