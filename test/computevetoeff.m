
clear 

reqAccVetoRate = 1.17e-05;
numHTrigs = 6443;
numXTrigs = 1697;
dbFileLoc = '/Users/pajith/working/dc/veto/omega_veto/test/sqlite/';

dbFileVec = {'L1-LSC-MICH_CTRL+L1-ASC-QPDX_P_L1-LSC-MICH_CTRL+L1-ASC-QPDX_P-data.db', ...   
    'L1-LSC-MICH_CTRL+L1-ASC-QPDX_Y_L1-LSC-MICH_CTRL+L1-ASC-QPDX_Y-data.db', ...
    'L1-LSC-PRC_CTRL+L1-ASC-QPDX_P_L1-LSC-PRC_CTRL+L1-ASC-QPDX_P-data.db', ...
    'L1-LSC-PRC_CTRL+L1-ASC-QPDX_Y_L1-LSC-PRC_CTRL+L1-ASC-QPDX_Y-data.db'};

dbFile = [dbFileLoc '/' dbFileVec{1}];

params = 'tau, r, rMax, trigHCentTime, trigXCentTime, trigHCentFreq, trigXCentFreq, trigHTrgSignf, trigXTrgSignf, trigHDuration, trigXDuration';

% get data from data base
mksqlite('open', dbFile);
Data = mksqlite(['select ' params ' from data']);
mksqlite('close');

tau     = [Data.tau];
r       = [Data.r];
rMax    = [Data.rMax];
trigHCentTime = [Data.trigHCentTime];
trigXCentTime = [Data.trigXCentTime];
trigHCentFreq = [Data.trigHCentFreq];
trigXCentFreq = [Data.trigXCentFreq];
trigHTrgSignf = [Data.trigHTrgSignf];
trigXTrgSignf = [Data.trigXTrgSignf];
trigHDuration = [Data.trigHDuration];
trigXDuration = [Data.trigXDuration];

timeShiftIdx = find(tau ~= 0);
zeroLagIdx   = find(tau == 0);

% figure
% plot(tau(timeShiftIdx), abs(r(timeShiftIdx)), 'b.')
% hold on
% plot(tau(zeroLagIdx), abs(r(zeroLagIdx)), 'rx')
% grid on


rTimeShift = abs(r(timeShiftIdx));
rZeroLag = abs(r(zeroLagIdx));
    
rVec = linspace(0,1,50);

[rVec, probDensTS] = calcprobdensity(rTimeShift, rVec);
[rVec, probDensZL] = calcprobdensity(rZeroLag, rVec);
                
numTimeShifts = length(unique(tau))-1;
analysisStartTime = min(trigHCentTime);
analysisEndTime = max(trigHCentTime);
    
% compute the veto threshold corresponding to the given accidental
% veto rate
[rThresh] = findvetothreshold(rTimeShift,  reqAccVetoRate, ...
                    numTimeShifts, analysisStartTime, analysisEndTime)
    
% recompute the accidental veto rate corresponding to the computed
% threshold -- sanity check
analysisDuration = analysisEndTime-analysisStartTime;
effNumSecs = analysisDuration*numTimeShifts;
N = length(find(rTimeShift >= rThresh));
accVetoRate = N/effNumSecs

% compute the veto efficiency, use percentage and accidental 
% veto prob.
vetoIdx = find(rZeroLag >= rThresh);
noVetoIdx = find(rZeroLag < rThresh);
vetoEff = 100*length(unique(trigHCentTime(zeroLagIdx(vetoIdx))))/numHTrigs;
usePerc = 100*length(unique(trigXCentTime(zeroLagIdx(vetoIdx))))/numXTrigs;
accVetoProb = length(find(rTimeShift >= rThresh))/(numHTrigs*numTimeShifts);

vetoEffCoinc = length(vetoIdx)/length(rZeroLag);

fprintf('# rThresh = %4.3f vetoEff = %4.3f accVetoProb = %3.2e accVetoRate = %3.2e\n', ...
    rThresh, vetoEff, accVetoProb, accVetoRate);

return

[vetoDensity] = hist(trigHCentTime(zeroLagIdx(vetoIdx)), timeVec);

figure
figure
stairs(timeVec-min(timeVec),vetoDensity, 'r')
grid on

deadTime = length(find(vetoDensity >= 1))


