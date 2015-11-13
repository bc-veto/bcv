function projUV = projectvector(u,v)
% 
% PROJECTVECTOR - find the projection of the (complex) vector v on to
% another vector u. 
%
% usage: projUV = projectvector(u,v)
%
% u, v   : vectors - real or complex. v and u must be of the same size
% projUV : projection of v on to u
% 
% P. Ajith, 15.03.05
% 
% $Id: projectvector.m,v 1.3 2006/09/21 14:05:24 ajith Exp $

% find the matrix size
matSize = size(v);
nCol    = matSize(2);

if size(u) ~= matSize
    error('### size of u and v should be the same');
end

% in the case of vectors
if min(matSize) == 1
    vDotU = sum(v.*conj(u));
    uSqr  = sum(u.*conj(u));
    projUV = (vDotU/uSqr).*u;
    
% in the case of matrices    
else
    vDotU = sum(v.*conj(u),1);
    uSqr  = sum(u.*conj(u),1);
    for iCol=1:nCol
        projUV(:,iCol) = (vDotU(iCol)./uSqr(iCol)).*u(:,(iCol));
    end
end
