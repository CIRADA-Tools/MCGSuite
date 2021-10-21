cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the definition of a particle
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module ParticleMod
      implicit none

      Type Particle
        real Flux               !The particle flux
        real Pos(0:2), Vel(0:2)                 !The position and velocity of each particle
        real AngPos(0:2)                    !The position of the particle in cylindrical coordinates
        real ProjectedPos(0:2), ProjectedVel(0:2)   !The position and velocity of the particles after projection
      end Type


      end module
