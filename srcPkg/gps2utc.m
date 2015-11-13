% GPS2UTC is a function that converts GPS seconds to UTC time.
%
% Usage: function UTC_time = gps2utc(GPS_time)
%
% GPS_time - Time in GPS seconds
%
% UTC_time = gps2utc(GPS_time)
% Converts GPS seconds to UTC time.
%
% Example:
%
% UTC_time = gps2utc(711129614)
% Will return: UTC_time = '2002-07-19 16:00:00'
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

function UTC_time = gps2utc(GPS_time)

% Set the GPS epoch.
GPS_Epoch = datenum('01-06-1980 00:00:00') * 86400;

% Calculate the numerical time from the GPS epoch and GPS time.
NUM_time = GPS_Epoch + GPS_time - 14;

% Calculate the UTC time and convert the result into a string.
UTC_time = datestr(NUM_time / 86400, 31);