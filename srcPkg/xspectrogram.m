function xspectrogram(tiling, transforms)
% XSPECTROGRAM Temporary function to display timefrequency spectrograms split over many Q planes
%
% XSPECTROGRAM displays a single time-frequency spectrogram of varying Q
% for XTILE time-frequency tiling.

% extract the number of planes
numberOfPlanes = tiling.numberOfPlanes;

% initialize the image
tf = [];

% loop over each Q plane
for plane = 1:numberOfPlanes,
    
    % ensure that there is one frequency per Q plane
    if(tiling.planes{plane}.numberOfRows ~= 1)
        error('more than one row in a plane; use WSPECTROGRAM');
    end
    
    % postpend the row to an image
    tf = [tf ; transforms{1}.planes{plane}.rows{1}.normalizedEnergies];
end

% display the image
imagesc(tf);
