#   Total number of galaxies in the sample
nTot=100

#   Masses
MassDict={'Flag':True,'FixedVal':9.5,'RanType':0,'Range':[7.5,8.5]}
#   Size Flag
BeamDict={'Flag':True,'FixedVal':5.,'RanType':1,'Range':[2.,8.]}
#   Inclination Flag
IncDict={'Flag':True,'FixedVal':45.,'RanType':0,'Range':[5.,90.]}
#   Postion Angles
PADict={'Flag':True,'FixedVal':0.,'RanType':0,'Range':[0.,360.]}
#   Velocity Dispersions
VelDispDict={'Flag':True,'FixedVal':8.,'RanType':0,'Range':[6.,14.]}
#   Decide if using the UDG mode
UDG_switch=True
VHIDict={'Flag':True,'FixedVal':23.0,'RanType':0,'Range':[15.,70.]}


n_Processors=10

# Output Folder
OutFolder='RandomSuite'
# Verbose plot switch
PlotVerbose=True
# Verbose file Switch
FileVerbose=True
#   Beta Configuration Options
#BetaConfigFile="Inputs.beta_config"

SuiteDict={'IncDict':IncDict,'MassDict':MassDict,'BeamDict':BeamDict,'PADict':PADict,'VelDispDict':VelDispDict,'nTot':nTot,'VHIDict':VHIDict}

