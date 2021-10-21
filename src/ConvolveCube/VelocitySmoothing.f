cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       fill a single ring with the particles necessary for
c       a tilted ring model
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module VelSmoothMod
      use BeamMod
      use DataCubeMod
      use CommonConsts


      implicit none

      contains

ccccccc
      subroutine SmoothVelChannels(DC,B)
      implicit none
      Type(DataCube), intent(INOUT) :: DC
      Type(Beam2D), INTENT(IN) :: B

      if(B%VelocitySmoothSwitch .eq. 0) then        !No velocity smoothing
        return
      elseif(B%VelocitySmoothSwitch .eq. 1) then    !Use Gaussian velocity smoothing
        call GaussSmoothCube(DC,B)          !In this file
      endif


      return
      end subroutine
cccccccc

ccccccccc
      subroutine GaussSmoothCube(DC,B)
      implicit none
      Type(DataCube), intent(INOUT) :: DC
      Type(Beam2D), INTENT(IN) :: B
      integer nChannelsToSmooth
      integer i, j, Indx(0:1)

c           Get the number of velocity channels to smooth
      nChannelsToSmooth=int(B%VelocitySmoothSigma/DC%DH%ChannelSize
     &                  *B%SigmaLengths)
c           Loop over all the pixels
      do i=0, DC%DH%nPixels(0)-1
        do j=0, DC%DH%nPixels(1)-1
c           Smooth the flux across all channels
            call GaussSmoothChannel(DC%DH
     &             ,DC%Flux(i,j,0:DC%DH%nChannels-1)
     &              ,nChannelsToSmooth,B%VelocitySmoothSigma)   !In this file
        enddo
      enddo

      return
      end subroutine
cccccccc

cccccccc
      subroutine GaussSmoothChannel(DH, FluxColumn
     &                  ,numSigmaChannels,Sigma)
      implicit none
      integer, INTENT(IN) :: numSigmaChannels
      real, INTENT(IN) :: Sigma
      Type(DataCubeHeader), INTENT(IN) :: DH
      real, INTENT(INOUT) :: FluxColumn(0:DH%nChannels-1)
c
      real TempColumn(0:DH%nChannels-1)
      integer i,j,k
      real deltaV,GaussFactor(-numSigmaChannels:numSigmaChannels)

      TempColumn=0.

      do i=-numSigmaChannels, numSigmaChannels
        deltaV=i*DH%ChannelSize
        GaussFactor(i)=1./sqrt((2.*Pi*Sigma**2.))
     &                  *exp(-(deltaV**2./(2*Sigma**2.)))
      enddo
      GaussFactor=GaussFactor/(sum(GaussFactor))

      do i=0, DH%nChannels-1
        do j=-numSigmaChannels,numSigmaChannels
            k=i+j
            if(k .ge. 0. .and. k .lt. DH%nChannels) then
                TempColumn(i)=TempColumn(i)
     &                  +GaussFactor(j)*FluxColumn(k)
            endif
        enddo
      enddo

      FluxColumn=TempColumn

      return
      end subroutine
ccccccccc

      end module
