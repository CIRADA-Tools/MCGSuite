cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines to calculate the moments 
c     of some a passed array of n-body particles per 
c     pixel element
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DataCubeMod
      implicit none

      Type DataCubeHeader
        integer nPixels(0:1),nChannels
        real PixelSize(0:1),ChannelSize
        real Start(0:2)
        character(8) AxisType(0:2),Units(0:2)
        integer DimensionUnitSwitch(0:3)    !0->PixelSize, 1->ChannelSize, 2->RefVal(0:1), 3->RefVal(2)
        integer UncertaintyUnitSwitch
        real Uncertainty
        character(10) FUnit,FType
        real Epoch
        integer PixelCenterIndx(0:1),ChannelCenterIndx
        real RefLocation(0:2), RefVal(0:2)
      end Type

      Type DataCube
        Type(DataCubeHeader) DH
        real,dimension(:,:),ALLOCATABLE :: Pixels
        real,dimension(:),ALLOCATABLE :: Channels
        real,dimension(:,:,:),ALLOCATABLE :: Flux

      end Type

      contains
ccccccc
      subroutine AllocateDataCube(DC)
c
      implicit none
      Type(DataCube),INTENT(INOUT) :: DC
      integer i, j,mPix
      real Delta

      mPix=maxval(DC%DH%nPixels)
      ALLOCATE(DC%Pixels(0:1,0:mPix-1))
      ALLOCATE(DC%Channels(0:DC%DH%nChannels-1))
      ALLOCATE(DC%Flux(0:DC%DH%nPixels(0)-1,0:DC%DH%nPixels(1)-1
     &          ,0:DC%DH%nChannels-1))
      DC%Flux=0.

      do j=0,1
        Delta=(0.0-(DC%DH%RefLocation(j)))
        DC%DH%Start(j)=DC%DH%RefVal(j)+Delta*DC%DH%PixelSize(j)

c        print*, "Allocating DC", j,DC%DH%Start(j),DC%DH%RefLocation(j)
c     &              ,DC%DH%RefVal(j)
c     &              ,DC%DH%nPixels(j),DC%DH%PixelSize(j)
        do i=0,mPix-1
            Delta=(real((i))-(DC%DH%RefLocation(j)))
            DC%Pixels(j,i)=DC%DH%RefVal(j)+Delta*DC%DH%PixelSize(j)
c            print*, "Pixel Check", j,i,DC%Pixels(j,i)
        enddo
      enddo

      j=2
      Delta=(0.0-(DC%DH%RefLocation(j)))
      DC%DH%Start(j)=DC%DH%RefVal(j)+(Delta
     &                  *DC%DH%ChannelSize)

c      print*, "Vel Start",DC%DH%Start(j),Delta
c     &          ,DC%DH%RefLocation(j),DC%DH%ChannelSize

      do i=0,DC%DH%nChannels-1
        Delta=(real((i))-(DC%DH%RefLocation(j)))
        DC%Channels(i)=DC%DH%RefVal(j)+Delta*DC%DH%ChannelSize
c        print*, i, DC%Channels(i), DC%DH%Start(2),DC%DH%ChannelSize
      enddo

c      DC%DH%PixelCenterIndx=DC%DH%nPixels/2
c      DC%DH%ChannelCenterIndx=DC%DH%nChannels/2
c      print*, "DataCube pixel (channel) center indices"
c     &          ,DC%DH%PixelCenterIndx
c     &          ,DC%DH%ChannelCenterIndx

      return
      end subroutine
cccccccc

cccccc
      subroutine DeAllocateDataCube(DC)
      implicit none
      Type(DataCube) DC

      DEALLOCATE(DC%Pixels)
      DEALLOCATE(DC%Channels)
      DEALLOCATE(DC%Flux)

      return
      end subroutine
cccccccc

      end module
