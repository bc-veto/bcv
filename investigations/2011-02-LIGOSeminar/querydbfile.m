function[MetaData, Data] = querydbfile(dbFile, params)

        mksqlite('open', dbFile);

        if (length(params)) >= 1
            Data = mksqlite(['select ' params ' from data']);
        else 
            Data = [];
        end
    
        query = 'select value from results where varName = ';
        numTrigsX = mksqlite([query, '''numTrigsX''']);
        numTrigsH = mksqlite([query, '''numTrigsH''']);
        SNRcutoffH = mksqlite([query, '''SNRcutoffH''']);
        SNRcutoffX = mksqlite([query, '''SNRcutoffX''']);
        couplingModel = mksqlite([query, '''couplingModel''']);
        rThresh = mksqlite([query, '''rThresh''']);
        highPassCutoff = mksqlite([query, '''highPassCutoff''']);
        analysisStartTime = mksqlite([query, '''analysisStartTime''']);
        analysisEndTime = mksqlite([query, '''analysisEndTime''']);
        analysisDuration = mksqlite([query, '''analysisDuration''']);
        reqAccVetoRate = mksqlite([query, '''reqAccVetoRate''']);
        safetyProbability = mksqlite([query, '''safetyProbability''']);
        vetoEfficiency = mksqlite([query, '''vetoEfficiency''']);
        accidentalVetoProb = mksqlite([query, '''accidentalVetoProb''']);
        efficiencyOverDeadtime = mksqlite([query, '''efficiencyOverDeadtime''']);
        deadTimePercentage = mksqlite([query, '''deadTimePercentage''']);
        totalInjectionNumber = mksqlite([query, '''totalInjectionNumber''']);
        vetoEfficiencyCoincTrigs = mksqlite([query, '''vetoEfficiencyCoincTrigs''']);
        accidentalVetoRate = mksqlite([query, '''accidentalVetoRate''']);
        vetoSignificance = mksqlite([query, '''vetoSignificance''']);
        Nvetoed = mksqlite([query, '''Nvetoed''']);
        deadTime = mksqlite([query, '''deadTime''']);
        logFile = mksqlite([query, '''logFile''']);
        configurationFile = mksqlite([query, '''configurationFile''']);
        mksqlite('close');
    
        MetaData.numTrigsX = numTrigsX.value;
        MetaData.numTrigsH = numTrigsH.value;
        MetaData.SNRcutoffH = SNRcutoffH.value;
        MetaData.SNRcutoffX = SNRcutoffX.value;
        MetaData.couplingModel = couplingModel.value;
        MetaData.rThresh = rThresh.value;
        MetaData.highPassCutoff = highPassCutoff.value;
        MetaData.analysisStartTime = analysisStartTime.value;
        MetaData.analysisEndTime = analysisEndTime.value;
        MetaData.analysisDuration = analysisDuration.value;
        MetaData.reqAccVetoRate = reqAccVetoRate.value;
        MetaData.safetyProbability = safetyProbability.value;
        MetaData.vetoEfficiency = vetoEfficiency.value;
        MetaData.accidentalVetoProb = accidentalVetoProb.value;
        MetaData.efficiencyOverDeadtime = efficiencyOverDeadtime.value;
        MetaData.deadTimePercentage = deadTimePercentage.value;
        MetaData.totalInjectionNumber = totalInjectionNumber.value;
        MetaData.vetoEfficiencyCoincTrigs = vetoEfficiencyCoincTrigs.value;
        MetaData.accidentalVetoRate = accidentalVetoRate.value;
        MetaData.vetoSignificance = vetoSignificance.value;
        MetaData.Nvetoed = Nvetoed.value;
        MetaData.deadTime = deadTime.value;
        MetaData.logFile = logFile.value;
        MetaData.configurationFile = configurationFile.value;

    


