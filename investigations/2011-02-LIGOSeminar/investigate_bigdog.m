
clear 

colTS = 'b';
colZL = 'r';

tL1 = 968654557.950;
tH1 = 968654557.957;

reqAccVetoRate = 1.65e-06;
numHTrigs = 1884;
numXTrigs = 1;

% location and file name of the sqlite database file
% dbFileLoc = 'db/H1';
% t0 = tH1;
% dbFileVec_H1;

dbFileLoc = 'db/L1';
t0 = tL1;
dbFileVec_L1;

% the parameters need to be retrieved from the database 
params = 'tau, r, rMax, trigHCentTime, trigXCentTime, trigHCentFreq, trigXCentFreq, trigHTrgSignf, trigXTrgSignf, trigHDuration, trigXDuration';

for iFile = 1:length(dbFileVec)

    % get data from data base
    dbFile = [dbFileLoc '/' dbFileVec{iFile}];
    mksqlite('open', dbFile);
    Data = mksqlite(['select ' params ' from data']);
    Data1 = mksqlite(['select value from results where varName = ''rThresh''']);
    mksqlite('close');
    
    % convert the resulting structures to be vectors 
    tau           = [Data.tau];
    r             = [Data.r];
    rMax          = [Data.rMax];
    trigHCentTime = [Data.trigHCentTime];
    trigXCentTime = [Data.trigXCentTime];
    trigHCentFreq = [Data.trigHCentFreq];
    trigXCentFreq = [Data.trigXCentFreq];
    trigHTrgSignf = [Data.trigHTrgSignf];
    trigXTrgSignf = [Data.trigXTrgSignf];
    trigHDuration = [Data.trigHDuration];
    trigXDuration = [Data.trigXDuration];
    rThresh       = [Data1.value];

    fprintf('%s: num_coinc_trigs = %d. ', dbFile, length(trigHCentTime));
    
    % find triggers near the time of the big dog. 
    bdIdx = find(trigHCentTime > t0-2 & trigHCentTime < t0+2);
    
    if length(bdIdx) > 0
    
        % identify triggers from the zero lag and time shift 
        timeShiftIdx = find(tau ~= 0);
        zeroLagIdx   = find(tau == 0);
        
        figure
        plot(tau(timeShiftIdx), abs(r(timeShiftIdx)), [colTS '.'])
        hold on
        plot(tau(zeroLagIdx), abs(r(zeroLagIdx)), [colZL '.'])
        plot(tau(bdIdx), abs(r(bdIdx)), 'gx')
        h1 = plot(tau, ones(size(tau))*rThresh, 'k-');
        set(h1, 'lineWidth', 1);
        grid on
        title(dbFile)
        
        rTimeShift = abs(r(timeShiftIdx));
        rZeroLag = abs(r(zeroLagIdx));
            
        rVec = linspace(0,1,50);
        
        [rVec, probDensTS] = calcprobdensity(rTimeShift, rVec);
        [rVec, probDensZL] = calcprobdensity(rZeroLag, rVec);
                        
        figure
        stairs(rVec, probDensZL, colZL)
        hold on
        stairs(rVec, probDensTS, colTS)
        h1 = plot([rThresh rThresh], [min(probDensZL) max(probDensZL)], 'k--');
        set(h1, 'lineWidth', 1);
        grid on
        xlabel('|r|')
        ylabel('Prob. density')
        title(dbFile)
    
    else
        fprintf('no trigger found near t0 \n');
    end
    
end
