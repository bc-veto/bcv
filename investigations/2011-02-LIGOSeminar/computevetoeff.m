function [vetoIdx, zeroLagIdx, rThresh, vetoEff, usePerc, accVetoProb, vetoSignific, vetoEffCoinc...
    ] = computevetoeff(tau, r, rMax, trigHCentTime, trigXCentTime, rThresh, chanName, ...
    reqAccVetoRate)
% 
% COMPUTEVETOEFF compute veto efficiency and other parameters using 
% the db file produced by the bilinear coupling veto code. 
% 
% usage: [vetoEff, usePerc, accVetoProb, vetoSignific, vetoEffCoinc, ...
%    tau, r, rMax] = computevetoeff(chanName, reqAccVetoRate)
%
% 
% P. Ajith, 6 October 2010
% 
% $Id:$

    colTS = 'b';
    colZL = 'r';

    numHTrigs = 1884;
    numXTrigs = 1;

    if length(r) > 0
    
        % identify triggers from the zero lag and time shift 
        timeShiftIdx = find(tau ~= 0);
        zeroLagIdx   = find(tau == 0);
        
        rTimeShift = abs(r(timeShiftIdx));
        rZeroLag = abs(r(zeroLagIdx));

        numTimeShifts = length(unique(tau))-1;
        analysisStartTime = min(trigHCentTime);
        analysisEndTime = max(trigHCentTime);
            
        % compute the veto threshold corresponding to the given accidental
        % veto rate
        [rThresh] = findvetothreshold(rTimeShift,  reqAccVetoRate, ...
                            numTimeShifts, analysisStartTime, analysisEndTime);
            
        figure
        subplot(121)
        plot(tau(timeShiftIdx), abs(r(timeShiftIdx)), [colTS '.'])
        hold on
        plot(tau(zeroLagIdx), abs(r(zeroLagIdx)), [colZL '.'])
        h1 = plot(tau, ones(size(tau))*rThresh, 'k--');
        set(h1, 'lineWidth', 1);
        grid on
        title(chanName)
        
            
        rVec = linspace(0,1,50);
        
        [rVec, probDensTS] = calcprobdensity(rTimeShift, rVec);
        [rVec, probDensZL] = calcprobdensity(rZeroLag, rVec);
                        
        subplot(122)
        stairs(rVec, probDensZL, colZL)
        hold on
        stairs(rVec, probDensTS, colTS)
        h1 = plot([rThresh rThresh], [min(probDensZL) max(probDensZL)], 'k--');
        set(h1, 'lineWidth', 1);
        grid on
        xlabel('|r|')
        ylabel('Prob. density')

        % recompute the accidental veto rate corresponding to the computed
        % threshold -- sanity check
        analysisDuration = analysisEndTime-analysisStartTime;
        effNumSecs = analysisDuration*numTimeShifts;
        N = length(find(rTimeShift >= rThresh));
        accVetoRate = N/effNumSecs;
        
        % compute the veto efficiency, use percentage and accidental 
        % veto prob.
        vetoIdx = find(rZeroLag >= rThresh);
        noVetoIdx = find(rZeroLag < rThresh);
        vetoEff = 100*length(unique(trigHCentTime(zeroLagIdx(vetoIdx))))/numHTrigs;
        usePerc = 100*length(unique(trigXCentTime(zeroLagIdx(vetoIdx))))/numXTrigs;
        accVetoProb = length(find(rTimeShift >= rThresh))/(numHTrigs*numTimeShifts);
        vetoSignific = vetoEff/100/accVetoProb;
        vetoEffCoinc = 100*length(vetoIdx)/length(rZeroLag);
        
        fprintf('rThresh = %2.1e vetoEff = %2.1e vetoEffCoinc = %2.1e\n# accVetoProb = %2.1e accVetoRate = %2.1e VetoSign = %2.1f\n', ...
            rThresh, vetoEff, vetoEffCoinc, accVetoProb, accVetoRate, vetoSignific);

    
    else
        fprintf('no trigger found \n');
    end
    
