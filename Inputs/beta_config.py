###         TILTED RING MODE PARAMETERS  -- will overwrite the default values
# Cloud mode
cmode = 0
# Base number density of particle clouds
CloudSurfDens = 1000
# random seed for the particle generation
ranseed=-5


####        PROFILES PARAMETERS  -- will overwrite the default values
nBinsPerRHI=100.
limR_RHI=1.


####        ADDITIONAL TILTED RING MODEL PARAMTERS
# x center (degrees)
RA=277.67543
# y center (degrees)
DEC=73.434828
# systemic velocity (km/s)
vsys=1403.93164636
# radial velocity (km/s)
vrad=0.
# derivative of the radial velocity as a function of height above zgradient height
dvdz=0.



###         DATA CUBE PARAMETERS
# Cube Shape  -- can overwrite the calculated cube shape
cube_shape=[200,200,200]
#   Beam flattening (relative to the FWHM defined in observatory_config.py)
Beam_Flattening=1.
#   Beam angle (in degrees)
Beam_PositionAngle=0.
# number of sigma lengths to use for beam smoothing -- overwrites default value of 5
nSigma=5.
# Switch for velocity smoothing (0=none, 1=Gaussian) [default is 0]
velocity_smooth_switch=0
# Width for the velocity smoothing if being used  (km/s) [default is channel width]
velocity_smooth_sigma=4.




