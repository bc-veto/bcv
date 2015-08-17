
clear 

colTS = 'b';
colZL = 'r';

reqAccVetoRate = 1.65e-06;
numHTrigs = 1884;
numXTrigs = 1;

% location and file name of the sqlite database file
% dbFileLoc = 'db/H1';
% t0 = tH1;
% dbFileVec_H1;

dbFileLoc = 'db/L1';
dbFileVec_L1;

% the parameters need to be retrieved from the database 
params = ['tau, r, rMax, trigHCentTime, trigXCentTime, trigHCentFreq,', ...
    'trigXCentFreq, trigHTrgSignf, trigXTrgSignf, trigHDuration, trigXDuration'];

%for iFile = 2:2:length(dbFileVec)
for iFile = 2:2:2

    % get data from data base
    dbFile1 = [dbFileLoc '/' dbFileVec{iFile}];
    dbFile2 = [dbFileLoc '/' dbFileVec{iFile-1}];

    %%%%%%%%%%%%%%%%%%%%%% pseudo channel - 1 %%%%%%%%%%%%%%%%%%%%%%%%%
    [MetaData1, Data1] = querydbfile(dbFile1, params);

    % convert the resulting structures to be vectors 
    tau1           = [Data1.tau];
    r1             = [Data1.r];
    rMax1          = [Data1.rMax];
    trigHCentTime1 = [Data1.trigHCentTime];
    trigXCentTime1 = [Data1.trigXCentTime];
    trigHCentFreq1 = [Data1.trigHCentFreq];
    trigXCentFreq1 = [Data1.trigXCentFreq];
    trigHTrgSignf1 = [Data1.trigHTrgSignf];
    trigXTrgSignf1 = [Data1.trigXTrgSignf];
    trigHDuration1 = [Data1.trigHDuration];
    trigXDuration1 = [Data1.trigXDuration];

    fprintf('\n#################################################### \n');
    fprintf('\n# channel - 1: %s\n', dbFile1, params);
    [vetoIdx1, zeroLagIdx1, rThresh1] = computevetoeff(tau1, r1, rMax1, trigHCentTime1,...
        trigXCentTime1, rThresh1, dbFile1, reqAccVetoRate);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%% pseudo channel - 2 %%%%%%%%%%%%%%%%%%%%%%%%%
    [MetaData2, Data2] = querydbfile(dbFile2, params);

    % convert the resulting structures to be vectors 
    tau2           = [Data2.tau];
    r2             = [Data2.r];
    rMax2          = [Data2.rMax];
    trigHCentTime2 = [Data2.trigHCentTime];
    trigXCentTime2 = [Data2.trigXCentTime];
    trigHCentFreq2 = [Data2.trigHCentFreq];
    trigXCentFreq2 = [Data2.trigXCentFreq];
    trigHTrgSignf2 = [Data2.trigHTrgSignf];
    trigXTrgSignf2 = [Data2.trigXTrgSignf];
    trigHDuration2 = [Data2.trigHDuration];
    trigXDuration2 = [Data2.trigXDuration];

    fprintf('\n# channel - 2: %s\n', dbFile2);
    [vetoIdx2, zeroLagIdx2, rThresh2] = computevetoeff(tau2, r2, rMax2, trigHCentTime2,...
        trigXCentTime2, rThresh2, dbFile2, reqAccVetoRate);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear idx r_comb1 r_comb2 r_comb3 rMax tau

%    %%%%%%%%%%%%%%%%%%%%%% combined channel %%%%%%%%%%%%%%%%%%%%%%%%%
    k = 1;
    for i=1:length(r1)
        idx = find(trigHCentTime2 == trigHCentTime1(i) & ...
            trigHCentFreq2 == trigHCentFreq1(i) & tau2 == tau1(i) & ...
            trigHTrgSignf2 == trigHTrgSignf1(i));

        if length(idx) >= 1
            x = r1(i);
            y = r2(idx(1));
            r_comb1(k) = x+1i*y;
            r_comb2(k) = x+y;
            r_comb3(k) = x*y;
            rMax(k) = rMax1(i) + 1i*rMax2(idx(1));
            tau(k) = tau1(i);
            k = k+1;
        end
    end

    fprintf('\n# combined r1 + i r2 \n')
    [vetoIdx, zeroLagIdx, rThresh_comb1] = computevetoeff(tau, abs(r_comb1), rMax, trigHCentTime2,...
        trigXCentTime2, rThresh2, [dbFile1 '-' dbFile2], reqAccVetoRate);

    fprintf('\n# combined r1 + r2 \n')
    [vetoIdx, zeroLagIdx] = computevetoeff(tau, abs(r_comb2), rMax, trigHCentTime2,...
        trigXCentTime2, rThresh2, [dbFile1 '-' dbFile2], reqAccVetoRate);

    fprintf('\n# combined r1 * r2 \n')
    [vetoIdx, zeroLagIdx] = computevetoeff(tau, abs(r_comb3), rMax, trigHCentTime2,...
        trigXCentTime2, rThresh2, [dbFile1 '-' dbFile2], reqAccVetoRate);

    r1Vec = linspace(min(abs(real(r_comb1))), max(abs(real(r_comb1))), 100);
    r2Vec = sqrt(rThresh_comb1^2 - r1Vec.^2);

    figure
    plot(abs(real(r_comb1(tau == 0))), abs(imag(r_comb1(tau == 0))), 'r.');
    grid on; hold on; 
    plot(abs(real(r_comb1(tau ~= 0))), abs(imag(r_comb1(tau ~= 0))), 'b.');
    h1 = plot([min(abs(real(r_comb1)))  max(abs(real(r_comb1)))], [rThresh2 rThresh2], 'g--')
    h2 = plot( [rThresh1 rThresh1], [min(abs(imag(r_comb1)))  max(abs(imag(r_comb1)))], 'g--')
    h3 = plot(r1Vec, r2Vec, 'k--')
    set(h1, 'linewidth', 1);
    set(h2, 'linewidth', 1);
    set(h3, 'linewidth', 1);
    xlabel('r_1', 'color', 'k')
    ylabel('r_2', 'color', 'k')

    r1Vec = linspace(min((real(r_comb1))), max((real(r_comb1))), 1000);
    r2Vec = sqrt(rThresh_comb1^2 - r1Vec.^2);

    figure
    plot((real(r_comb1(tau == 0))), (imag(r_comb1(tau == 0))), 'r.');
    grid on; hold on; 
    plot((real(r_comb1(tau ~= 0))), (imag(r_comb1(tau ~= 0))), 'b.');
    h1 = plot([min(real(r_comb1))  max(real(r_comb1))], [rThresh2 rThresh2], 'g--', ...
        [min(real(r_comb1))  max(real(r_comb1))], [-rThresh2 -rThresh2], 'g--')
    h2 = plot( [rThresh1 rThresh1], [min(imag(r_comb1))  max(imag(r_comb1))], 'g--', ...
        [-rThresh1 -rThresh1], [min(imag(r_comb1))  max(imag(r_comb1))], 'g--')
    h3 = plot(r1Vec, r2Vec, 'k--', r1Vec, -r2Vec, 'k--')
    set(h1, 'linewidth', 1);
    set(h2, 'linewidth', 1);
    set(h3, 'linewidth', 1);
    xlabel('r_1', 'color', 'k')
    ylabel('r_2', 'color', 'k')



%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
