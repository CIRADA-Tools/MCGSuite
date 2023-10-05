cccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module contains the routines for writing out a
c       data cube to a file.
c
ccccccccccccccccccccccccccccccccccccccccccccccccc

      module DataCubeOutputsMod
      use DataCubeMod
      use BeamMod


      implicit none

      contains

ccccccc
c           This routine is the main routine for writing out a FITS cubes
c               Note that many of these routines are adapted from the cfitsio cookbook
      subroutine WriteDataCubeToFITS(DC,Beam,fname,Name)
      implicit none
      Type(DataCube), INTENT(IN):: DC
      Type(Beam2D),INTENT(IN) :: Beam
      character(*), INTENT(IN) :: fname
      character(*),INTENT(IN):: Name

      integer status,unit
      integer naxes(3)
      print*, "Outputting Data Cube"

      call InitializeFITSFileIO(fname,status,unit)
      call WriteFITSHeader(DC,Beam,fname,status,unit,naxes
     &          ,Name)
      call WriteFITSBody(DC,fname,status,unit,naxes)


c       Finally close the file using the cfitsio routines
C  The FITS file must always be closed before exiting the program.
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)


      return
      end subroutine
cccccccc




ccccc
      subroutine InitializeFITSFileIO(fname,status,unit)
      implicit none
      character(*), INTENT(IN) :: fname
      integer, INTENT(OUT) :: status,unit
      integer blocksize
c           Set the status parameter for the FITSIO routines
      status=0
      call ftgiou(unit,status)      !Get the unit for the IO
c       Check if the file already exists and delete it
      call deletefile(fname,status)
C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      call ftinit(unit,fname,blocksize,status)
c      print*, "Initialized Fits unit", unit, trim(fname)

      return
      end subroutine
cccccc

ccccc
c           This routine writes out all the necessary header information
      subroutine WriteFITSHeader(DC,Beam,fname,status,unit,naxes
     &              ,Name)
      use CommonConsts
      implicit none
      Type(DataCube), INTENT(IN):: DC
      Type(Beam2D),INTENT(IN) :: Beam
      character(*), INTENT(IN) :: fname
      character(*),INTENT(IN) :: Name
      integer, INTENT(IN) :: unit
      integer, INTENT(INOUT) :: status
      integer, INTENT(INOUT) :: naxes(3)

      logical simple, extend
      integer naxis,bitpix

      character(10) date,time,zone
      character(20) tString
      integer date_I(8)
      integer test1,test2

      integer naxesT(3),naxisT

      simple=.true.
      naxis=3
      extend=.true.

      bitpix=-32
      naxes(1:2)=DC%DH%nPixels(0:1)
      naxes(3)=DC%DH%nChannels

      print*, "Outputting Data Cube", naxes

      naxesT=naxes
      naxisT=naxis

      test1=0
      test2=1
c      print*, "Consistency header check", unit, simple
c     &      ,bitpix,naxis,naxes,extend,status
C  Write the required header keywords to the file
      call ftphpr(unit,simple,bitpix,naxisT
     &      ,naxesT,test1,test2,extend,status)
c      call ftphps(unit,bitpix,naxis
c     &      ,naxes,status)
c      print*, "Writing header", naxis,naxes


c       Write the 'optional' keywords that are also necessary for the file
      call ftpkys(unit,'CTYPE1',DC%DH%AxisType(0),
     &          'PRIMARY AXIS Type ',status)
      call ftpkye(unit,'CRPIX1',DC%DH%RefLocation(0)+1
     &          ,8,
     &          'X-ref-pixel',status)
      call ftpkye(unit,'CDELT1',DC%DH%PixelSize(0)/3600.,8,
     &          'X-increment',status)
      call ftpkye(unit,'CRVAL1',
     &          DC%DH%RefVal(0)/3600.,8,
     &          'Value at X-ref-pixel',status)
      call ftpkys(unit,'CUNIT1',DC%DH%Units(0),
     &          ' ',status)


      call ftpkys(unit,'CTYPE2',DC%DH%AxisType(1),
     &          'PRIMARY AXIS Type ',status)
      call ftpkye(unit,'CRPIX2',DC%DH%RefLocation(1)+1
     &          ,8,
     &          'Y-ref-pixel',status)
      call ftpkye(unit,'CDELT2',DC%DH%PixelSize(1)/3600.,8,
     &          'Y-increment',status)
      call ftpkye(unit,'CRVAL2',
     &          DC%DH%RefVal(1)/3600.,8,
     &          'Value at Y-ref-pixel',status)
      call ftpkys(unit,'CUNIT2',DC%DH%Units(0),
     &          ' ',status)

      call ftpkys(unit,'CTYPE3',DC%DH%AxisType(2),
     &          'PRIMARY AXIS Type ',status)
      call ftpkye(unit,'CRPIX3',DC%DH%RefLocation(2)+1
     &          ,8,
     &          'Z-ref-channel',status)
      call ftpkye(unit,'CDELT3',DC%DH%ChannelSize*1000.,8,
     &          'Z-increment',status)
      call ftpkye(unit,'CRVAL3',
     &          DC%DH%RefVal(2)*1000.,8,
     &          'Value at Z-ref-channel',status)
      call ftpkys(unit,'CUNIT3',DC%DH%Units(2),
     &          ' ',status)


      call ftpkye(unit,'EPOCH',DC%DH%Epoch,8,
     &          'EPOCH ',status)
      call ftpkys(unit,'BUNIT',DC%DH%FUnit,
     &          'no comment ',status)
      call ftpkys(unit,'BTYPE',DC%DH%FType,
     &          ' ',status)

c       The beam needs to be converted to degrees from pixel
c           units
      call ftpkye(unit,'BMAJ',
     &          Beam%BeamMajorAxis
     &          *DC%DH%PixelSize(1)/3600,8,
     &          'Major Axis of the Beam',status)

      call ftpkye(unit,'BMIN',
     &          Beam%BeamMinorAxis
     &          *DC%DH%PixelSize(1)/3600,8,
     &          'Major Axis of the Beam',status)

      call ftpkye(unit,'BPA',
     &          Beam%BeamPositionAngle*180./Pi,8,
     &          'Major Axis of the Beam',status)



      call Date_and_Time(date,time,zone,date_I)
      write(tString,'(I2.2,"-",I2.2, "-", I4.4)')
     &          date_I(3),date_I(2),date_I(1)

      call ftpkys(unit,'DATE',tString,
     &          ' ',status)

      call ftpkys(unit,'ORIGIN',"MockCubeGenerator",
     &          ' ',status)

c      call ftpkys(unit,'OBJECT',Name,
c     &          ' ',status)

c       Header values needed for CARTA?
      call ftpkye(unit,'RestFreq',
     &           1.42040575179E+09,8,
     &          'Rest frequency of HI in Hz',status)
      call ftpkys(unit,'SPECSYS',"BARYCENT",
     &          ' ',status)

      call ftpkys(unit,'RADESYS',"FK5",
     &          ' ',status)

c       Add a bunch of extra bits to deal with WCS definitions

      call ftpkye(unit,'PC01_01',1.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC02_01',0.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC03_01',0.0,8,
     &          ' ',status)

      call ftpkye(unit,'PC01_02',0.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC02_02',1.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC03_02',0.0,8,
     &          ' ',status)

      call ftpkye(unit,'PC01_03',0.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC02_03',0.0,8,
     &          ' ',status)
      call ftpkye(unit,'PC03_03',1.0,8,
     &          ' ',status)

      call ftpkye(unit,'LONPOLE',180.0,8,
     &          ' ',status)
      call ftpkye(unit,'LATPOLE',0.0,8,
     &          ' ',status)

      return
      end subroutine
ccccccc




cccccc

      subroutine WriteFITSBody(DC,fname,status,unit,naxes)
      implicit none
      Type(DataCube), INTENT(IN):: DC
      character(*), INTENT(IN) :: fname
      integer, INTENT(IN) :: unit,naxes(3)
      integer, INTENT(INOUT) :: status

      integer*8 group,fpixel,nelements
      real,ALLOCATABLE:: FTest(:,:,:)
      integer i,j
    

c Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument; in this case the 'J' indicates that the array has an
C  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.

c      print*, "Writing Fits Body"
      ALLOCATE(FTest(naxes(1),naxes(2),naxes(3)))
      do i=1, naxes(1)
        do j=1,naxes(2)
            FTest(i,j,1:naxes(3))=DC%Flux(i-1,j-1,0:naxes(3)-1)
        enddo
      enddo
    

      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
c      print*, "check on shape/size", nelements, group,fpixel
c     &          ,shape(FTest)
c      call ftppre(unit,group,fpixel,nelements,DC%Flux,status)
      call ftppre(unit,group,fpixel,nelements,FTest,status)
c      call ffppxe(unit,group,fpixel,nelements,FTest,status)
c      print*, "Done Fits body write"

      DEALLOCATE(FTest)

      return
      end subroutine
cccccc



C *************************************************************************
c           The necessity of this routine was indicated in the cfitsio library.
c               Rather than checking on the existence of the file and deleting if
c               necessary, I've re-written it to just delete the file if it already
c               exists as we don't need to append anything to the file
      subroutine deletefile(filename,status)
      implicit none
C  A simple little routine to delete a FITS file
      character(*), INTENT(IN) :: filename
      integer, INTENT(INOUT) :: status
      character(500) Script
      logical FileExists

      inquire(file=trim(filename),exist=FileExists)
      if(FileExists) then
        Script="rm "//trim(filename)
c        print*, Script
        call system(Script)
      endif



      return
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end module
