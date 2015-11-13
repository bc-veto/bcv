% UTC2GPS is a function that converts UTC time to GPS seconds.
%
% Usage: function GPS_time = utc2gps(UTC_time)
%
% UTC_time - Coordinated Universal Time
%
% GPS_time = utc2gps(UTC_time)
% Converts UTC time to GPS seconds.
% UTC_time can also be an array of UTC times.
%
% Examples:
%
% GPS_time = utc2gps('2002-07-19 16:00:00')
% Will return: GPS_time = 711129614
%
% GPS_time = utc2gps(['2002-07-19 16:00:00'; '2001-07-19 16:00:00'])
% Will return: GPS_time = [711129614 ; 679593614]
%
% M. Hewitson
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 07-07-09
% P. Ajith, 30-06-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Vetoes for Binary-Coalescence Triggers Using Known Instrumental     %
%                                Couplings                                %
%                                                                         %
%                      Modified By: Aaron B. Pearlman                     %
%                        Mentor: Ajith Parameswaran                       %
%                        Creation Date: 22-06-2009                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GPS_time = utc2gps(UTC_time)

% Set the GPS epoch.
GPS_Epoch = datenum('01-06-1980 00:00:00') * 86400;

% Get the number of entries in UTC_time.
numberOfEntries = size(UTC_time, 1);

% Allocate space for the array of GPS times to be returned.
GPS_time = zeros(1, numberOfEntries);

% Loop over the number of entries in UTC_time.
for entry = 1 : numberOfEntries
    % Get each UTC time from the UTC_time array.
	CurrUTC = UTC_time(entry, :);
    
	% Reformt string to MATLAB format MM-DD-YYY.
	CurrUTC = strcat(CurrUTC(6 : 10), '-', CurrUTC(1 : 4), CurrUTC(11 : length(CurrUTC)));
    
    % Calculate the numerical time from the selected UTC time.
	NUM_time = datenum(CurrUTC) * 86400;
    
    % Convert the result into GPS seconds using the GPS epoch.
	GPS_time(entry) = round(NUM_time - GPS_Epoch + 14);
end