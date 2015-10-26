import numpy as np

class TrigStruct:
  def __init__(self, trigStartTime, trigEndTime , trigCentralTime, trigCentFreq,
	       trigSign):
    self.startTime = trigStartTime
    self.endTime = trigEndTime
    self.centralTime = trigCentralTime
    self.centralFrequency = trigCentFreq
    self.triggerSignificance = trigSign
