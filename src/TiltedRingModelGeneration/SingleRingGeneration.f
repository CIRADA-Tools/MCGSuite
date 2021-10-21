cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module SingleRingGenerationMod
      use ParticleMod
      use TiltedRingMod


      implicit none

      contains

ccccccc
      subroutine Ring_CalcNumParticles(R,pixelarea,cmode,CloudSurfDens)
      use CommonConsts
      implicit none
      Type(Ring), INTENT(INOUT):: R
      real, INTENT(IN) :: pixelarea   !The pixel area is needed to get the
      integer,INTENT(IN) ::cmode
      real, INTENT(IN) :: CloudSurfDens
      real Pixel_Ring
      real DensMultiplications
      real Rl,Rh

c           Get the rough area of the ring in pixels
c      Pixel_Ring=2.*Pi*R%Rmid*R%Rwidth /pixelarea      !Original in 3DBarolo GalMod
      Rl=R%Rmid-R%Rwidth/2.
      Rh=R%Rmid+R%Rwidth/2.
      Pixel_Ring=Pi*(Rh**2.-Rl**2)!/pixelarea
c       Calculate a term to get the rough number of clouds per pixel area.  It is normalized by
c           the surface density so that each particle has roughly the same amount of flux
      DensMultiplications=CloudSurfDens
     &                      *((R%Sigma)**real(cmode))
c           Get an integer number of particles
      R%nParticles=int(DensMultiplications*Pixel_Ring)+1
c      print*, "Single Ring Position Angle",R%PositionAngle*180./Pi
c      print*, "Single Ring Center",R%CentPos
c      print*, "Number of particles", R%nParticles,R%Sigma
c      print*, R%Rwidth,R%Rmid,R%Sigma,CloudSurfDens,cmode
c     &          ,Pixel_Ring,DensMultiplications, R%PositionAngle*180./Pi


      return
      end subroutine
cccccccc

ccccccc
      subroutine Ring_ParticleGeneration(R,idum)
      use BasicRanNumGen
      use CommonConsts
      implicit none
      Type(Ring), INTENT(INOUT) :: R     !The ring that we are generating particles for
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      integer i
      real Rmin,Rmax, Area       !The minimum and maximum radii of the ring and the ring area

c           To randomly sample the ring area we need Rmin and Rmax
      Rmin=R%Rmid-R%Rwidth/2.       !We need Rmin and Rmax for the ring first
      Rmax=R%Rmid+R%Rwidth/2.
      Area=Pi*(Rmax**2.-Rmin**2)
c      print*, "Single Ring Area Check", Area,Rmin,Rmax

c      print*, "MRing Params", R%CentPos,R%Inclination, R%PositionAngle
c     &          ,R%VSys, R%VRot, R%VRad, R%VDisp,R%Vvert
c     &          ,R%dvdz,R%Sigma,R%z0,R%zGradiantStart
c      print*, "Consistency Check", R%VSys,R%VRad
c           Now go through all the particles
      do i=0, R%nParticles-1
        call Ring_ParticlePosSelect(R,idum, Rmin,Rmax,i)
        call ParticlePosProject(R%P(i), R%Inclination, R%PositionAngle)
        call ParticlePos_NewCenter(R%P(i), R%CentPos)
        call Ring_CalcParticle_VSys(R,i,idum)
        call Ring_CalcParticleFlux_Basic(R,i, Area)

c        print*, "Random Check", R%P(i)%AngPos, R%P(i)%Pos
c        print*, "Projected Random Check", R%P(i)%ProjectedPos(0:1)
c        print*, "PRojected Velocity", R%P(i)%ProjectedVel(2)
      enddo

c      print*, "particle Flux ChecK", sum(R%P(0:R%nParticles-1)%Flux)


      return
      end subroutine
cccccccc


ccccccc
      subroutine Ring_ParticlePosSelect(R,idum, Rmin, Rmax, PartID )
      use BasicRanNumGen
      use CommonConsts
      implicit none
      Type(Ring), INTENT(INOUT) :: R     !The ring that we are generating particles for
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      integer, INTENT(IN) :: PartID
      real, INTENT(IN) :: Rmin,Rmax
      real RR, Theta, Z     !Temporary cylindrical coordinates of each particle

c           First get a random radius -- however this is not uniform in R.  We want
c           equal area sampling we'll be using a sqrt distribution.
      RR=ran2(idum)*(Rmax**2.-Rmin**2.)             !Sample a uniform squared radius in the correct range
      RR=sqrt(RR+Rmin**2.)          !Add on the minimum radius squared and take the sqrt
c           Next get a random angle
      Theta=ran2(idum)*2.*Pi
c           Finally get a height from a sech^2 distribution (use atanh on a random distrubtion and scale
c                   by the height
      Z=atanh((2.*ran2(idum)-1.))*R%z0

c           Store the cylindrical coordinates
      R%P(PartID)%AngPos(0)=RR
      R%P(PartID)%AngPos(1)=Theta
      R%P(PartID)%AngPos(2)= Z
c           Store the Cartesian Coordinates
      R%P(PartID)%Pos(0)=RR*cos(Theta)
      R%P(PartID)%Pos(1)=RR*sin(Theta)
      R%P(PartID)%Pos(2)= Z

      return
      end subroutine
ccccccccc

cccccccc
      subroutine ParticlePosProject(P,Inclination,PositionAngle)
      implicit none
      Type(Particle), INTENT(INOUT) :: P
      real, INTENT(IN) :: Inclination, PositionAngle
      real cpa,spa,XTemp,YTemp
c       First get the inclined positions in the major axis frame (aligned with X-axis)
c           in the temporary arrays
      XTemp=P%Pos(0)
      YTemp=P%Pos(1)*cos(Inclination)+P%Pos(2)*sin(Inclination)
c       Now rotate the positions by the position angle
      cpa=cos(PositionAngle)        !Get the cos and sin values
      spa=sin(PositionAngle)
c       Calculate the coordinates in the projected array
      P%ProjectedPos(0)=XTemp*cpa-YTemp*spa
      P%ProjectedPos(1)=XTemp*spa+YTemp*cpa

c      print*, P%ProjectedPos, PositionAngle*180./3.14, XTemp,YTemp
      return
      end subroutine
cccccccccccc

cccccccc
      subroutine ParticlePos_NewCenter(P,NewCent)
      implicit none
      Type(Particle), INTENT(INOUT) :: P
      real, INTENT(IN) :: NewCent(0:1)

c       Shift the projected position to lie about the new center
c      print*, "Initial Pos", P%ProjectedPos,NewCent
      P%ProjectedPos(0)=P%ProjectedPos(0)+NewCent(0)
      P%ProjectedPos(1)=P%ProjectedPos(1)+NewCent(1)
c      print*, "Final Pos", P%ProjectedPos
      return
      end subroutine
cccccccccccc

cccccc
      subroutine Ring_CalcParticle_VSys(R,PartID,idum)
      use BasicRanNumGen
      implicit none
      Type(Ring), INTENT(INOUT) :: R
      integer, INTENT(IN) :: PartID
      integer, INTENT(INOUT) :: idum          !idum is a seed for the random number generator
      real VRotP,VRadP,VVertP,VDispP
      real V_FromRotation, V_FromRadial,V_FromVertical
      real delR,delZ
      real cTheta,sTheta
      integer DerivativeSide

c           First get the particle's rotational, radial, and vertical motions
c           at the current radial positions

c               To do the linear interpolation, figure out which side you're on
c      delR=(R%P(PartID)%AngPos(0)-R%Rmid)
c      if(delR .le. 0) then
c        DerivativeSide=0
c      else
c        DerivativeSide=1
c      endif

c      delZ=abs(R%P(PartID)%AngPos(2))-R%zGradiantStart
      VRotP=R%VRot !+ R%dVRot_dR(DerivativeSide)*delR
      VRadP=R%VRad !+ R%dVRad_dR(DerivativeSide)*delR
      VDispP=R%VDisp !+ R%dVDisp_dR(DerivativeSide)*delR
c           Check if the particle is far enough above the plane that V_rot is affected
      if(abs(R%P(PartID)%AngPos(2)) .gt. R%zGradiantStart) then
        VRotP=VRotP-R%dvdz*delZ
      endif
c
c           Now get the sin and cos of the angle
c
      cTheta=cos(R%P(PartID)%AngPos(1))
      sTheta=sin(R%P(PartID)%AngPos(1))
c
c           Next get the observed velocity due to the main components
c
      V_FromRotation=VRotP*cTheta*sin(R%Inclination)
      V_FromRadial=VRadP*sTheta*sin(R%Inclination)
      V_FromVertical=R%Vvert*cos(R%Inclination)         !Note the Vvert is the global vertical motion of a ring
c        Sum up the various components
      R%P(PartID)%ProjectedVel(2)=R%VSys+V_FromRotation
     &                      +V_FromRadial+V_FromVertical
c
c           Finally get a random velocity from the dispersion and add it in to the projected velocity
c
      R%P(PartID)%ProjectedVel(2)=R%P(PartID)%ProjectedVel(2)
     &                      +gasdev(idum)*VDispP

c      print*, "hmmm", VDispP, R%VDisp

      return
      end subroutine
ccccccccc

cccccc
      subroutine Ring_CalcParticleFlux_Basic(R,PartID,Area)
      implicit none
      Type(Ring), INTENT(INOUT) :: R
      integer, INTENT(IN) :: PartID
      real,INTENT(IN):: Area
      real SigmaP, del R
      integer DerivativeSide
c           Get the local Sigma using the derivative and radius
c      delR=(R%P(PartID)%AngPos(0)-R%Rmid)
c      if(delR .le. 0) then
c        DerivativeSide=0
c      else
c        DerivativeSide=1
c      endif

c      SigmaP=R%Sigma+R%dSigma_dR(DerivativeSide)*delR
      SigmaP=R%Sigma

      R%P(PartID)%Flux=SigmaP*Area/real(R%nParticles)

c      print*, "Flux calc check", SigmaP,Area,
c     &                  R%P(PartID)%Flux,R%nParticles

      return
      end subroutine
cccccccc

      end module
