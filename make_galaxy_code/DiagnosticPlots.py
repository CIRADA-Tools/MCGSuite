#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from . import StandardRelations as SR
from . import ObjectDefinitions as OD
from . import CubeAnalysis as CA

import gc

def MakeAllPlots(GalaxyIO,Galaxy,DataCube,Profiles,TiltedRing):
    if GalaxyIO.plot_verbose == False:
        print("Not making any plots")
        return
    
    BasePlotParams={'font.size': 15,'axes.linewidth':2
            ,'xtick.major.size':6,'xtick.minor.size':3
            ,'xtick.major.width':1,'xtick.minor.width':1
            ,'ytick.major.size':6,'ytick.minor.size':3
            ,'ytick.major.pad':10,'xtick.major.pad':10
            ,'ytick.major.width':1,'ytick.minor.width':1
            ,'xtick.labelsize':17 ,'ytick.labelsize':17
            ,'axes.labelsize': 25}

    mpl.rcParams.update(BasePlotParams)


    RHI_AS=SR.KpcToArcSec(Galaxy.RHI,Galaxy.distance)
    MakeMomentMapsPlot(GalaxyIO,Galaxy.pa,DataCube.beam_fwhm,DataCube.noise,RHI_AS,Galaxy.inclination,Galaxy.distance,Galaxy.VHI,TiltedRing.vsys,Galaxy.veldisp,GalaxyIO.MapPlotName)

    MakeProfilesPlots(Profiles,Galaxy,GalaxyIO.ProfilePlotName)

    return

def MakeProfilesPlots(Profiles,Galaxy,PlotName):
    label_size = 25.0
    num_size = 21.0
    dashthick = 1.5
    axlim = 1.1
    Rlim = axlim*Profiles.limR_RHI*Galaxy.RHI
    Vlim = axlim*Galaxy.VHI
    SBlim = 12.0
    fig, ((axtop),(axbot)) = plt.subplots(2,1,figsize=(10,15))


    axtop.set_xlim(0,Rlim)
    axtop.set_ylim(0,Vlim)
    axtop.tick_params(axis='both', which='major', labelsize=num_size)
    axtop.set_xlabel('R (kpc)',fontsize=label_size)
    axtop.set_ylabel('V (km/s)',fontsize=label_size)
    axtop.text(1.02*Galaxy.RHI,0.05*Galaxy.VHI,"$\mathrm{R_{HI}}$",fontsize=label_size, c='tab:orange')
    axtop.text(1.15*Galaxy.RHI,0.92*Galaxy.VHI,"$\mathrm{V_{HI}}$",fontsize=label_size, c='tab:blue')
    axtop.plot([0.,Rlim],[Galaxy.VHI,Galaxy.VHI],'--',linewidth=dashthick)
    axtop.plot([Galaxy.RHI,Galaxy.RHI],[0.,Vlim],'--',linewidth=dashthick)
    axtop.plot(Profiles.R,Profiles.Vrot,linewidth=2.5)
    
    axbot.set_xlim(0,Rlim)
    axbot.set_ylim(0.3,SBlim)
    axbot.tick_params(axis='both', which='major', labelsize=num_size)
    axbot.set_xlabel('R (kpc)',fontsize=label_size)
    axbot.text(1.02*Galaxy.RHI,0.35,"$\mathrm{R_{HI}}$",fontsize=label_size, c='tab:orange')
    axbot.text(1.15*Galaxy.RHI,0.75,"$\mathrm{N_{HI}}$",fontsize=label_size, c='tab:blue')
    axbot.set_ylabel('$\Sigma \, \mathrm{(M_\odot / pc^2)}$',fontsize=label_size)
    axbot.set_yscale("log")
    axbot.plot([0.,Rlim],[1.0,1.0],'--',linewidth=dashthick)
    axbot.plot([Galaxy.RHI,Galaxy.RHI],[0.,SBlim],'--',linewidth=dashthick)
    axbot.plot(Profiles.R,Profiles.SB,linewidth=2.5)

    plt.savefig(PlotName,bbox_inches='tight')

    plt.close(fig)

    gc.collect()

def MakeMomentMapsPlot(GalaxyIO,PositionAngle,BeamFWHM,noise,RHI,Inc,distance,VHI,VSys,VDisp,PlotName):
    #GetProfilesAndMaps(GalaxyIO)
    PX,PY,Vels,VelocityProfile,Mom0,Mom1,Mom2,SourceProfile,SourceMom0,SourceMom1,SourceMom2,PVMajor,PVMinor,PVMajorSource,PVMinorSource=CA.GetProfilesAndMaps(GalaxyIO,PositionAngle,BeamFWHM,noise,distance)
    
    
    PXX,PYY=np.meshgrid((PX-PX[int(len(PX)/2)])*3600.,(PY-PY[int(len(PY)/2)])*3600.)
    PXX2,VS=np.meshgrid((PY-PY[int(len(PY)/2)])*3600.,Vels)
    #print((PY-PY[int(len(PY)/2)])*3600.)

    fig = plt.figure(figsize=(20,10))
    
    left=0.15
    base=0.15
    w=0.3
    h=w*2
    buf=0.12
    wbuf=0.12
    hbuf=0.2
    
    cMap='plasma'
    MomentLabels=(r'RA offset (")','DEC offset (")')
    PVLabels=(r'offset (")','v (km s$^{-1}$)')
    VelLabels=(r'v (km s$^{-1}$)','F')
    
    ax1=fig.add_axes([left,base+2*(h+hbuf),w,h])
    ax2=fig.add_axes([left+(w+wbuf),base+2*(h+hbuf),w,h])
    ax3=fig.add_axes([left+2*(w+wbuf),base+2*(h+hbuf),w,h])
    
    ax4=fig.add_axes([left,base+h+hbuf,w,h])
    ax5=fig.add_axes([left+(w+wbuf),base+h+hbuf,w,h])
    ax6=fig.add_axes([left+2*(w+wbuf),base+h+hbuf,w,h])
    
    #ax7=fig.add_axes([left,base,w,h])
    #ax8=fig.add_axes([left+(w+wbuf),base,w,h])
    #ax9=fig.add_axes([left+2*(w+wbuf),base,w,h])


#ax10=fig.add_axes([left,base-(h+hbuf),w,h])
#ax11=fig.add_axes([left+(w+wbuf),base-(h+hbuf),w,h])
#ax12=fig.add_axes([left+2*(w+wbuf),base-(h+hbuf),w,h])
    
    delP=(PX[1]-PX[0])*3600
    CentVal=[0,0]
    SLims=[0.,1.1*np.max(Mom0)]
    VLims=[VSys-1.1*VHI,VSys+1.1*VHI]
    DispersionPlotMax=np.max([1.1*VDisp,1.1*np.nanmax(Mom2)])   #To deal with projection effects on the Mom2 map
    VDLims=[0.,DispersionPlotMax]
    
    
    
    
    MomentMap_Ellipse(ax1,Mom0,cMap,PXX,PYY,MomentLabels,'Mom0',2*RHI,Inc,PositionAngle,CentVal,SLims)
    MomentMap_Ellipse(ax2,Mom1,cMap,PXX,PYY,MomentLabels,'Mom1',2*RHI,Inc,PositionAngle,CentVal,VLims)
    MomentMap_Ellipse(ax3,Mom2,cMap,PXX,PYY,MomentLabels,'Mom2',2*RHI,Inc,PositionAngle,CentVal,VDLims)

 
 #MomentMap_Ellipse(ax7,SourceMom0,cMap,PXX,PYY,MomentLabels,'Source Mom0',2*RHI,Inc,PositionAngle,CentVal,SLims)
 #MomentMap_Ellipse(ax8,SourceMom1,cMap,PXX,PYY,MomentLabels,'Source Mom1',2*RHI,Inc,PositionAngle,CentVal,VLims)
 #MomentMap_Ellipse(ax9,SourceMom2,cMap,PXX,PYY,MomentLabels,'Mom2',2*RHI,Inc,PositionAngle,CentVal,VDLims)
    
    
    MomentMap(ax4,PVMajor,cMap,PXX2,VS,PVLabels,'PV-Major')
    MomentMap(ax5,PVMinor,cMap,PXX2,VS,PVLabels,'PV-Minor')
    VelocityProfilePlot(ax6,Vels,VelocityProfile,VelLabels,'Profile')
  
  
  #MomentMap(ax10,PVMajorSource,cMap,PXX2,VS,PVLabels,'Source PV-Major')
  #MomentMap(ax11,PVMinorSource,cMap,PXX2,VS,PVLabels,'Source PV-Minor')
  #VelocityProfilePlot(ax12,Vels,SourceProfile,VelLabels,'Source Profile')

    plt.savefig(PlotName,bbox_inches='tight')

    plt.close(fig)

    del fig, Mom0, Mom1, Mom2, PVMajor,PVMinor,VelocityProfile
    del SourceMom0,SourceMom1,SourceMom2,PVMajorSource,PVMinorSource,SourceProfile

    gc.collect()

def MomentMap(ax,Map,cMap,X,Y,Labels,Title):
    #print("Moment map shapes", np.shape(Map), np.shape(X),np.shape(Y))
    ax.pcolormesh(X,Y,Map,cmap=cMap,shading='auto')
    FormatPlot(ax,Labels,Title)

def MomentMap_Ellipse(ax,Map,cMap,X,Y,Labels,Title,a,inc,pa,center,CLims):
    ax.pcolormesh(X,Y,Map,cmap=cMap,vmin=CLims[0],vmax=CLims[1],shading='auto')
    ellip=np.cos(inc*np.pi/180.)
    b=a*ellip
    angle=(90+pa)
    Ell=Ellipse(center, a, b, angle,edgecolor='cyan',facecolor='none',lw=5)
    ax.add_patch(Ell)
    FormatPlot(ax,Labels,Title)



def VelocityProfilePlot(ax,Vels,Flux,Labels,Title):
    ax.plot(Vels,Flux,ls='-',color='k')
    FormatPlot(ax,Labels,Title)

def FormatPlot(ax,Labels,Title):
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.set_xlabel(Labels[0])
    ax.set_ylabel(Labels[1])
    ax.text(0.5, 1.05,Title,fontsize=35,horizontalalignment='center',transform=ax.transAxes)

