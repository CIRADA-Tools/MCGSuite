cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains routines for simple rotations using the
c       particle class
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module RotationsMod
      use ParticleMod

      implicit none

      contains

cccccccc
c       Rotate the position and velocity of some particle using a given
c           rotation matrix
c
      subroutine ParticleRotate(P, rows,cols, RotationMatrix)
      implicit none
      integer rows,cols
      Type(Particle) P
      real RotationMatrix(0:rows-1,0:cols-1)        !Indexing being done to match C
      integer i

c       Rotate the positions and store in the projected position vector
      call RotateVector(P%Pos,RotationMatrix,P%ProjectedPos)
c         Rotate the velocities and store in the projected velocity vector
      call RotateVector(P%Vel,RotationMatrix,P%ProjectedVel)

      return
      end subroutine
ccccccc

cccccc
c           Rotate a 3-element vector using a rotation matrix and store it in
c           rotated vector
c
      subroutine RotateVector(Vector,RotationMatrix,RotatedVector)
      implicit none
      real Vector(0:2), RotatedVector(0:2)      !Index using C standard
      real RotationMatrix(0:2,0:2)              !Index using C standard
      integer i,j
c
      RotatedVector(0:2)=0.     !Set the rotated vector to zero initially
      do i=0,2                  !  Loop through columns
        do j=0,2                !  Loop through rows
            RotatedVector(j)=RotatedVector(j)
     &                          +RotationMatrix(j,i)*Vector(i)      !Do the matrix multiplication
        enddo
      enddo

      return
      end subroutine
cccccccc


ccccccc
c           Calculate a 3x3 rotation matrix corresponding to an inclination and
c           position angle
c
      subroutine CalculateRotationMatrix(Inclination, PositionAngle
     &                  ,rows,cols,RotationMatrix)
      implicit none
      real Inclination, PositionAngle
      integer rows,cols
      real RotationMatrix(0:rows-1,0:cols-1)        !Indexing to match C

      real ci,si,cpa,spa         !The cos and sin for the inclination and position angle

c      Calculate the sin and cos of both angles
      ci=cos(Inclination)
      si=sin(Inclination)
      cpa=cos(PositionAngle)
      spa=sin(PositionAngle)

c         Use a XZY Euler matrix for the rotations into the projected plane
c       The X rotation angle=inclination, Z rotation angle=position angle, Y rotation angle=0
c       Note that X-Y is spatial and v_z will be the radial velocity
      RotationMatrix(0,0)= cpa;
      RotationMatrix(0,1)= -spa;
      RotationMatrix(0,2)= 0.;

      RotationMatrix(1,0)=ci*spa ;
      RotationMatrix(1,1)= ci*cpa;
      RotationMatrix(1,2)= -si;

      RotationMatrix(2,0)= si*spa ;
      RotationMatrix(2,1)= si*cpa;
      RotationMatrix(2,2)= ci;

      return
      end subroutine
cccccccc

      end module
