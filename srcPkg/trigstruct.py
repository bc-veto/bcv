import numpy as np
# Structure which holds the information about the triggers

class TrigStruct:
  def __init__(self, trigStartTime, trigEndTime , trigCentralTime, trigCentFreq,
	       trigSign):
    self.startTime = trigStartTime
    self.endTime = trigEndTime
    self.centralTime = trigCentralTime
    self.centralFrequency = trigCentFreq
    self.triggerSignificance = trigSign
