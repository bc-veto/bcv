function [vetoEff, usePerc, accVetoProb, accVetoRate, vetoEffCoinc, rThresh] = applyvetothreshold(resultsDir,... 
    reqAccVetoRate, gpsTriggerHList, numXTrigs)

            rVec = linspace(0,1,50);

            resultsFile = [resultsDir,'/','corrstat_timeshift.dat']

            % if there is an existing file with this name, rename it 
            if fopen(resultsFile) >= 0
                system(['mv ', resultsDir, '/corrstat_timeshift.dat ', ...
                    resultsDir, '/corrstat_timeshift.old']);
			end

            % cat data from different time shifts to a single file
            system(['cat ', resultsDir, '/','corrstat_timeshift*.dat > ', resultsFile]);

            if fopen(resultsFile) >= 0

              dat = load(resultsFile);

              if ~isempty(dat)
                k = 1;
                tau  = dat(:,k); k = k+1;
                r    = dat(:,k); k = k+1;
                rMax = dat(:,k); k = k+1;
              	trigHCentTime = dat(:,k); k = k+1;
                trigXCentTime = dat(:,k); k = k+1;
                trigHCentFreq = dat(:,k); k = k+1;
                trigXCentFreq = dat(:,k); k = k+1;
                trigHTrgSignf = dat(:,k); k = k+1;
                trigXTrgSignf = dat(:,k); k = k+1;
                trigHDuration = dat(:,k); k = k+1;
                trigXDuration = dat(:,k); k = k+1;
    
                % plot the histograms of the time-shift and zero lag
                timeShiftIdx = find(tau ~= 0);
                zeroLagIdx   = find(tau == 0);
                
                rTimeShift = abs(r(timeShiftIdx));
                rZeroLag = abs(r(zeroLagIdx));
    
                [rVec, probDensTS] = calcprobdensity(rTimeShift, rVec);
                [rVec, probDensZL] = calcprobdensity(rZeroLag, rVec);
                
                numTimeShifts = length(unique(tau));
                analysisStartTime = min(trigHCentTime);
                analysisEndTime = max(trigHCentTime);
    
                % compute the veto threshold corresponding to the given accidental
                % veto rate
                [rThresh] = findvetothreshold(rTimeShift,  reqAccVetoRate, ...
                    numTimeShifts, analysisStartTime, analysisEndTime);
    
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
                numHTrigs = length(gpsTriggerHList.centralTime);
                vetoEff = length(unique(trigHCentTime(zeroLagIdx(vetoIdx))))/numHTrigs;
                usePerc = length(unique(trigXCentTime(zeroLagIdx(vetoIdx))))/numXTrigs;
                accVetoProb = length(find(rTimeShift >= rThresh))/length(rTimeShift);
                vetoEffCoinc = length(vetoIdx)/length(rZeroLag);

				
        else % no coincident triggers
          vetoEff = nan; 
          usePerc = nan;
          accVetoProb = nan;
          accVetoRate = nan;
          vetoEffCoinc = nan;
          rThresh = nan;
        end
    else
      error('### unable to open the results file.');
    end
