% DEBUGTIMESERIES - Create time-series plots of the raw and high-passed
% data from channels H and X, along with an identification of the central
% times of the coincident triggers.
%
% Aaron B. Pearlman <aaronp1@umbc.edu>, 06-08-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Vetoes for Binary-Coalescence Triggers Using Known Instrumental     %
%                                Couplings                                %
%                                                                         %
%                      Modified By: Aaron B. Pearlman                     %
%                        Mentor: Ajith Parameswaran                       %
%                        Creation Date: 06-08-2009                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
% 1. Plot the raw data - Channel H.

plotString = 'H1MICHCTRL-';
chanHNameString = strrep(chanHName, '_', '\_');

trigIdxH = find(rawTimeChH >= trigHCentTime);

% Plot the raw data and shift to time 0.
figure(1)
subplot(2, 1, 1)
plot(rawTimeChH - segStartTime, rawDataChH, 'b')
hold on
plot(rawTimeChH(trigIdxH(1)) - segStartTime, rawDataChH(trigIdxH(1)), 'ko')
plot(rawTimeChH(1) - segStartTime, 0, 'k*')
plot(rawTimeChH(end) - segStartTime, 0, 'k*')
grid on
title(sprintf('%s - Raw Data Channel H', chanHNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Data')

%**************************************************************************
% 2. Plot the high-passed data - Channel H.

% Plot the high-passed data and shift to time 0.
figure(1)
subplot(2, 1, 2)
plot(rawTimeChH - segStartTime, rawDataChHHP, 'r')
hold on
plot(rawTimeChH(trigIdxH(1)) - segStartTime, rawDataChHHP(trigIdxH(1)), ...
    'ko')
plot(rawTimeChH(1) - segStartTime, 0, 'k*')
plot(rawTimeChH(end) - segStartTime, 0, 'k*')
grid on
title(sprintf('%s - High-Passed Data Channel H', chanHNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Data')

saveas(gcf, sprintf('%sTimeSeriesChH.fig', plotString), 'fig')
saveas(gcf, sprintf('%sTimeSeriesChH.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-TimeSeriesChH.ps'

%**************************************************************************
% 3. Plot the raw data - Channel X.

chanXNameString = strrep(chanXName, '_', '\_');

trigIdxX = find(rawTimeChX >= trigXCentTime + timeShift);

% Plot the raw data and shift to time 0.
figure(2)
subplot(2, 1, 1)
plot(rawTimeChX - segStartTime, rawDataChX, 'b')
hold on
plot(rawTimeChX(trigIdxX(1)) - segStartTime, ...
    rawDataChX(trigIdxX(1)), 'ko')
plot(rawTimeChX(1) - segStartTime, 0, 'k*')
plot(rawTimeChX(end) - segStartTime, 0, 'k*')
grid on
title(sprintf('%s - Raw Data Channel X', chanXNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Data')

%**************************************************************************
% 4. Plot the high-passed data - Channel X.

% Plot the high-passed data and shift to time 0.
figure(2)
subplot(2, 1, 2)
plot(rawTimeChX - segStartTime, rawDataChXHP, 'r')
hold on
plot(rawTimeChX(trigIdxX(1)) - segStartTime, ...
    rawDataChXHP(trigIdxX(1)), 'ko')
plot(rawTimeChX(1) - segStartTime, 0, 'k*')
plot(rawTimeChX(end) - segStartTime, 0, 'k*')
grid on
title(sprintf('%s - High-Passed Data Channel X', chanXNameString))
xlabel(['Time Since ', num2str(segStartTime), ' [s]'])
ylabel('Data')

saveas(gcf, sprintf('%sTimeSeriesChX.fig', plotString), 'fig')
saveas(gcf, sprintf('%sTimeSeriesChX.pdf', plotString), 'pdf')
print -dpsc 'H1MICHCTRL-TimeSeriesChX.ps'

%**************************************************************************