   !!======================================================================
   !!                       ***  MODULE  icedrv_init_SIMBA ***
   !!=====================================================================
   !! History :  1.0   !  10-18  created by M. Plante 
   !! 			            based on code in CONCEPTS
   !!----------------------------------------------------------------------  
   !! Gets the forcing data from GEM, at a specific location, to force Icepack
   !!
   !!           If we initialise with buoy data, the location is taken 
   !!	        from the buoy data. Otherwise, the location is hardcoded
   !!           in icedrv_init.F90, to a location close to Nain.
   !!
   !!----------------------------------------------------------------------
   !!   prepare_gem_forcing     : Compute the atmospheric fluxes over ice
   !!----------------------------------------------------------------------!=========================================================================
!
! Get the initial conditions use to initialize Icepack
! from SIMBA data
!        
!=========================================================================
!       History : 1.0 ,  Nov. 2018, author: Mathieu Plante 
!

MODULE icedrv_init_SIMBA
  
    use icedrv_kinds 
    
    IMPLICIT NONE  
    PUBLIC 
    
    integer :: n, ns, ni, nsdata, nidata
    character(char_len_long), public :: & 
         data_buoy_dir           ! top directory for forcing data
    character*512 filename
    real, public, allocatable, dimension(:) ::  Tsnow, Tice    
    real (kind=dbl_kind), public, save :: lat_buoy, lon_buoy    
CONTAINS     

    SUBROUTINE get_buoy_data(nslyr, nilyr, Tair_buoy, hs, hi, Tns, Tni)

      implicit none

      INTEGER, INTENT(IN) :: nslyr, nilyr
 
      real (kind=dbl_kind) , intent(out) :: Tair_buoy, hs, hi
      real (kind=dbl_kind) :: lat_col, lon_col
      real (kind=dbl_kind), dimension(nslyr), intent(out)  ::  Tns 
      real (kind=dbl_kind), dimension(nilyr), intent(out)  ::  Tni     
      real  ::  zns, zni, x1, x2, f1, f2      
      real(kind=dbl_kind), allocatable :: profile(:) 

      ! print *, 'OPENING FILE?' 
      
      ! filename = data_buoy_dir !'/home/map005/data/eccc-ppp2/SIMBA_data/gca0103td2017022715'
      filename = data_buoy_dir !'/home/mwang/ICEPACK_MP/Forcing/SIMBA_2024_init.txt'

      print *, filename
         OPEN(UNIT=1,FILE=filename,FORM="FORMATTED",STATUS="OLD",ACTION="READ")
      print *, 'reading is not the problem...'
         ! read(1,FMT=*) lat_col
	 ! print *, 'nor is reading the latitude...'       
         ! read(1,FMT=*) lon_col
	 ! print *, lat_buoy, lon_buoy

         ! read(1,FMT=*) Tair_buoy
         ! read(1,FMT=*) hs
         ! read(1,FMT=*) hi
         ! read(1,FMT=*) nsdata 
         ! read(1,FMT=*) nidata  
        
         lat_col = 54.96069
         lon_col = -59.67123
         hs = 30 ! cm
         hi = 39 ! cm
         Tair_buoy = -15 ! deg. C
         nsdata = 10 ! Snow temperature profile (sensors 80-89)
         nidata = 19 ! Ice temperature profile (sensors 90-108)

         hs = hs/1d2
         hi = hi/1d2
         lat_buoy = lat_col
	 lon_buoy = lon_col+360.0
         print *, 'lon : ', lon_buoy, 'lat : ', lat_buoy
         print *, 'Tair_buoy : ', Tair_buoy, 'hs : ', hs, 'hi : ', hi      
         print *, 'allocating the profiles' 
         
         ALLOCATE( Tsnow(nsdata))                 
         ALLOCATE( Tice(nidata) )  !  , STAT=ierror 
         
         ! -----------------------------------------------------
         ! Read the full temperature profile from the file (all 241 sensors)
         ! Added my MW, Mar. 7
        print *, 'Allocating full temperature profile for 241 sensors'
        ALLOCATE(profile(241))
        
        ! Reading the profile data from file
        do n = 1, 241
            read(1,*) profile(n)    ! Read temperatures from file
        enddo

        ! Extract snow data (indices 80 to 89 from the profile)
        Tsnow = profile(80:89)

        ! Extract ice data (indices 90 to 108 from the profile)
        Tice = profile(90:108)

        ! Print extracted profiles for validation
        print *, 'Snow temperature profile (80-89): ', Tsnow
        print *, 'Ice temperature profile (90-108): ', Tice

        ! ------------------------------------------------------

         
        do ns = 1, nslyr
	    zns = ((2d0*ns)-1)*(hs*1d2)/(2d0*nslyr)
	    x1 = floor(zns/2d0)
	    x2 = x1+1
	    f1 = Tsnow(x1)
	    f2 = Tsnow(x2)	
	    Tns(ns) = f1 + ((zns/2d0 - x1)/(x2-x1))*(f2-f1)
	    print *, 'interpolation of internal temp : ', zns, x1, x2, f1, f2, Tns(ns)
	enddo
	
        do ni = 1, nilyr
	    zni = ((2d0*ni)-1)*(hi*1d2)/(2d0*nilyr)
	    x1 = floor(zni/2d0)
	    x2 = x1+1
	    f1 = Tice(x1)
	    f2 = Tice(x2)	
	    Tni(ni) = f1 + ((zni/2d0 - x1)/(x2-x1))*(f2-f1)   
	    print *, 'interpolation of internal temp : ', zni, x1, x2, f1, f2, Tni(ni)
	enddo
	
	print *, 'Tns = ', Tns
        print *, 'Tni = ', Tni
      return	
      END SUBROUTINE get_buoy_data

    END MODULE icedrv_init_SIMBA
