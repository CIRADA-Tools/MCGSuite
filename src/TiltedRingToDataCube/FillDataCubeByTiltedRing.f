cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines needed to
c       build a tilted ring model full of particles.
c          The routines assume the tilted ring has been allocated
c           and the individual ring parameters have been set.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module FillDataCubeWithTiltedRingMod
      use TiltedRingMod
      use DataCubeMod
      use ParticleMod


      implicit none
      contains

ccccc
      subroutine FillDataCubeWithTiltedRing(DC,TR)
      implicit none
      Type(DataCube), INTENT(INOUT) :: DC
      Type(TiltedRingModel), INTENT(IN) :: TR

      integer i,j,k,CellIndex(0:2)
      logical InBounds

      Type(Particle) PTest

      print*, "Filling DataCube"
c           Make sure the initial cube is empty
c      print*, "Cube center", DC%DH%Start(0:2)
      DC%Flux=0.
c           Loop through all particles
      do i=0,TR%nRings-1
        do j=0, TR%R(i)%nParticles-1
c                   Get the index for the cell where the particle lives
            call FindParticleCellLocation(TR%R(i)%P(j),DC,CellIndex)
c                   Check that the particle is inside the cube
            call CheckIfInCube(CellIndex,DC,InBounds)
c                   If it is in the cube, add the flux to the proper cell
            if(InBounds) then
                DC%Flux(CellIndex(0),CellIndex(1),CellIndex(2))=
     &                  DC%Flux(CellIndex(0),CellIndex(1),CellIndex(2))+
     &                  TR%R(i)%P(j)%Flux
            endif
        enddo
      enddo


c      print*, "Filled data cube flux check", sum(DC%Flux)
c      print*, "1st part flux", TR%R(0)%P(0)%Flux

c      print*, "TR Center Position Check"
c      PTest%ProjectedPos(0:1)=TR%R(0)%CentPos(0:1)
c      PTest%ProjectedVel(2)=TR%R(0)%VSys-4.
c      print*, "Center Position",PTest%ProjectedPos(0:1)
c     &          ,PTest%ProjectedVel(2)
c      call FindParticleCellLocation(PTest,DC,CellIndex)
c      print*, "Cell Index Calc", CellIndex
c      print*, "X Position Check", PTest%ProjectedPos(0)
c     &          , DC%DH%Start(0), DC%DH%PixelSize(0)
c     &          ,PTest%ProjectedPos(0)-DC%DH%Start(0)
c     & ,(PTest%ProjectedPos(0)-DC%DH%Start(0))/DC%DH%PixelSize(0)
c     &          ,PTest%ProjectedPos(0)/3600., DC%DH%Start(0)/3600.

c      print*, "Y Position Check", PTest%ProjectedPos(1)
c     &          , DC%DH%Start(1), DC%DH%PixelSize(1)
c     &          ,PTest%ProjectedPos(1)-DC%DH%Start(1)
c     & ,(PTest%ProjectedPos(1)-DC%DH%Start(1))/DC%DH%PixelSize(1)
c     &          ,PTest%ProjectedPos(1)/3600., DC%DH%Start(1)/3600.

      return
      end subroutine
cccccc

cccccc
      subroutine FindParticleCellLocation(P,DC,CellIndex)
      implicit none
      Type(Particle), INTENT(IN) :: P
      Type(DataCube), INTENT(IN) :: DC
      integer, INTENT(OUT) :: CellIndex(0:2)
      integer i,j

      do i=0,1
        CellIndex(i)=int(P%ProjectedPos(i)+0.5)
c-DC%DH%Start(i))
c     &                  /(DC%DH%PixelSize(i)) )

c        print*, "Cell Indx Check", i,CellIndex(i)
c     &          ,P%ProjectedPos(i)
c     &          ,DC%DH%PixelSize(i)
      enddo
      CellIndex(2)=int((P%ProjectedVel(2)-DC%DH%Start(2))
     &                  /(DC%DH%ChannelSize) +0.5)

      i=2


      return
      end subroutine
ccccccc

ccccccc
      subroutine CheckIfInCube(CellIndex,DC,InBounds)
      implicit none
      Type(DataCube), INTENT(IN) :: DC
      integer, INTENT(IN) :: CellIndex(0:2)
      logical, INTENT(OUT) :: InBounds
      integer i

      InBounds=.True.
c           Check the spatial dimensions
      do i=0,1
        if(CellIndex(i) .lt. 0 .or.
     &                  CellIndex(i) .ge. DC%DH%nPixels(i)) then
            InBounds=.False.
        endif
      enddo
c           Check the velocity
      if(CellIndex(2) .lt. 0 .or.
     &                  CellIndex(2) .ge. DC%DH%nChannels) then
        InBounds=.False.
      endif
c      print*, "Inbounds", InBounds, CellIndex
c     &              ,DC%DH%nPixels, DC%DH%nChannels


      return
      end subroutine
cccccccc


      end module
