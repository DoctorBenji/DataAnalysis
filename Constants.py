# In vivo constants
RecoveryTimepoints = {'PRE': 1,'IMM': 2, 1: 3, 2: 4, 5: 5, 10: 6}
HzFrequencies = [1, 5, 8, 10, 15, 20, 40, 50, 80, 100, 150, 200]
ISO_Percents = [10, 20, 30, 40, 50, 60, 70, 80]
END_ROM = -18.99


# Single fibre constants
pCas = ['4.5', '5.5', '5.7', '6.2', '6.4', '6.6', '7.0']


# Background info for test:
SampleRate:int = 10000
Stiffness_Time:int = 0.9 * SampleRate
ktr_start_time:int = int(1.016 * SampleRate)
ktr_end_time:int = int(6 * SampleRate)
PowerLoads = [10, 20, 30, 40, 50, 60, 70, 80]
ForceClampTimes = [(i * SampleRate, i * SampleRate + 4999) for i in (19.5, 20.0, 20.5, 21.0)]

Torque = r'Torque (mN$\cdot$m)'

