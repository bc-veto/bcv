
for 

    mksqlite('open', dbFile);

    Data(i).numTrigsX = mksqlite(['select value from results where varName = ''numTrigsX''']);
    Data(i).numTrigsH = mksqlite(['select value from results where varName = ''numTrigsH''']);
    Data(i).SNRcutoffH = mksqlite(['select value from results where varName = ''SNRcutoffH''']);
    Data(i).SNRcutoffX = mksqlite(['select value from results where varName = ''SNRcutoffX''']);
    Data(i).couplingModel = mksqlite(['select value from results where varName = ''couplingModel''']);
    Data(i).rThresh = mksqlite(['select value from results where varName = ''rThresh''']);
    Data(i).highPassCutoff = mksqlite(['select value from results where varName = ''highPassCutoff''']);
    Data(i).analysisStartTime = mksqlite(['select value from results where varName = ''analysisStartTime''']);
    Data(i).analysisEndTime = mksqlite(['select value from results where varName = ''analysisEndTime''']);
    Data(i).analysisDuration = mksqlite(['select value from results where varName = ''analysisDuration''']);
    Data(i).reqAccVetoRate = mksqlite(['select value from results where varName = ''reqAccVetoRate''']);
    Data(i).safetyProbability = mksqlite(['select value from results where varName = ''safetyProbability''']);
    Data(i).vetoEfficiency = mksqlite(['select value from results where varName = ''vetoEfficiency''']);
    Data(i).accidentalVetoProb = mksqlite(['select value from results where varName = ''accidentalVetoProb''']);
    Data(i).efficiencyOverDeadtime = mksqlite(['select value from results where varName = ''efficiencyOverDeadtime''']);
    Data(i).deadTimePercentage = mksqlite(['select value from results where varName = ''deadTimePercentage''']);
    Data(i).totalInjectionNumber = mksqlite(['select value from results where varName = ''totalInjectionNumber''']);
    Data(i).vetoEfficiencyCoincTrigs = mksqlite(['select value from results where varName = ''vetoEfficiencyCoincTrigs''']);
    Data(i).accidentalVetoRate = mksqlite(['select value from results where varName = ''accidentalVetoRate''']);
    Data(i).vetoSignificance = mksqlite(['select value from results where varName = ''vetoSignificance''']);
    Data(i).Nvetoed = mksqlite(['select value from results where varName = ''Nvetoed''']);
    Data(i).deadTime = mksqlite(['select value from results where varName = ''deadTime''']);
    Data(i).logFile = mksqlite(['select value from results where varName = ''logFile''']);
    Data(i).configurationFile = mksqlite(['select value from results where varName = ''configurationFile''']);

    mksqlite('close');
