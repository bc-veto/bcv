function [rThresh] = findvetothreshold(rTimeShift,  reqAccVetoRate, ...
    numTimeShifts, analysisStartTime, analysisEndTime)
% 
% FINDVETOTHRESHOLD - find the threshold on the cross-correlation statistic
% r corresponding to a required accidental veto rate. 
% 
% usage: [rThresh] = findvetothreshold(rTimeShift,  reqAccVetoRate, ...
%                       numTimeShifts, analysisStartTime, analysisEndTime)
% 
% rTimeShift        - a vector of r values from the time-shift
% reqAccVetoRate    - required accidental veto rate
% numTimeShifts     - number of timeshifts employed
% analysisStartTime - start time of the ananlysis (secs)
% analysisEndTime   - end time of the analysis (secs)
% rThresh           - threshold on r corresponding to the required acc. veto rate
% 
% P. Ajith, 14 Aug 2009
% 
% $Id: findvetothreshold.m 162 2009-08-14 20:11:29Z ajith $

if analysisStartTime > analysisEndTime
    error('analysisStartTime is larger than analysisEndTime');
end

% effective number of seconds in the data
effNumSecs = (analysisEndTime-analysisStartTime)*numTimeShifts;

Options = optimset('MaxIter', 1e6,...
                   'MaxFunEvals', 1e6,...
                   'TolFun', 1e-122,...
                   'Display', 'on',...
                   'Diagnostics', 'on');

rThreshInit = 0.5;
lowBound = -1.;
uppBound = 1.;

rThresh = fzero(@(rThresh ) accvetoratediff(rThresh, rTimeShift, effNumSecs, reqAccVetoRate), ...
    rThreshInit,  Options);

%threshVec = linspace(-1,1,1e3);
%for i=1:length(threshVec)
%    [deltaVec(i)] = accvetoratediff(threshVec(i), rTimeShift, effNumSecs, reqAccVetoRate);
%end
%
%figure
%plot(threshVec,deltaVec,'r')
%grid on


% compute the accidental veto rate corresponding to a threshold and subtract it
% from the required accidental veto rate. 
function [delta] = accvetoratediff(rThresh, rTimeShift, effNumSecs, reqAccVetoRate)

    rThresh;        

    reqAccVetoRate;

    % number of triggers with r > rThresh
    N = length(find(rTimeShift >= rThresh));

    % compute the accidental veto rate corresponding to the given threshold 
    accVetoRate = N/effNumSecs;

    % find the difference btwn the computed acc veto rate and the required rate 
    delta = accVetoRate-reqAccVetoRate;
