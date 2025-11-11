MODULE gemdrv_read_rpn
   !!======================================================================
   !!         ***  from CONCEPTS MODULE  fldread_rpn  ***
   !! Read input field for surface boundary condition 
   !!                 (rpn files)
   !!=====================================================================
   !! History CONCEPTS:  9.0  !  08-07  (F. Roy with original subroutines 
   !!                                           provided by  M. Desgagne) 
   !!  	      ICEPACK:   0.0  !  11-2018  (M. Plante, with subroutines from above) 
   !!----------------------------------------------------------------------
   !!   fld_read_rpn : 
   !!                   read input fields used for the computation of the
   !!                   surface boundary condition (written in rpn files)
   !!----------------------------------------------------------------------

    USE icedrv_kinds  
    USE icedrv_calendar
    USE gemdrv_sbc_rpntls

   IMPLICIT NONE
   PUBLIC   

   !! * Routine accessibility
   
   INTEGER :: iun_std(1000), &  ! Unit number vector
  &           nfile_std, &      ! Number of std files
  &           kt_sbc            ! Time step used for atm. data interpolation, or treatment of coupling 
   LOGICAL :: atm_file_L

   REAL(wp) ::   dt_gem_atm !: spacing between forcing fields (GEM) (hours at this stage)
   REAL(wp) ::   dt_gem_prc !: precip dt used for cumulation in GEM (hours at this stage)
                            !: equal to dt_gem_atm if no average in pre-process
                    
   LOGICAL  ::   ln_gem_intrp = .TRUE.   !: switch to activate linear time-interpolation
   INTEGER  ::   nh_gem_offs = 0          !: offset to match dates when searching forcing fields
                    !: example, if daily averages correspond to hour 3, nh_gem_offs = 3
                    !: must be positive
   LOGICAL  ::   ln_gem_avgrad = .FALSE.   !: switch to activate 24h averaging for radiation
   INTEGER, PUBLIC, DIMENSION(2) :: lev_gem_forc=(/ 28257976, 11950 /) 
                    !: Vertical level descriptor (equivalent to ip1), the first
                    !: one is the newstyle version, and the second one the
                    !: oldstyle. It allows to mix oldstyle and newstyle code
                    !: in atmospheric forcing std files.
                    !: The oldstyle code (second value) will only be used
                    !: if fields coded in newstyle are not found.
                    !: The default here is eta=0.995 
   TYPE, PUBLIC ::   FLD        !: Input field related variables
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) ::   fnow       ! input fields interpolated to now time step
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:,:) ::   fdta       ! 2 consecutive record of input fields
   END TYPE FLD  
   
   REAL(wp), PUBLIC, SAVE ::   &
                 zlev_gem   = 2       !: optional fixed forcing level from namelist
                                        !: allow to skip reading of PX,P0
                                        !: do nothing if left to -9.
  

   INTEGER, PUBLIC  ::   jpi = 1801   ! = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   !: first  dimension
   INTEGER, PUBLIC  ::   jpj =  1251 ! = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   !: second dimension  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: lat_rpn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: lon_rpn   
   
   REAL(wp), PUBLIC, PARAMETER ::   RGASV  =.46151e+3      ! J K-1 kg-1; gas constant 
                                                           ! for water vapour
   REAL(wp), PUBLIC, PARAMETER ::   CAPPA  =.28549121795   ! RGASD/CPD         ! ; Von Karman constant
   REAL(wp), PUBLIC, PARAMETER ::   RGASD  =.28705e+3    
   REAL(wp), PUBLIC, PARAMETER ::   DELTA  =.6077686814144 ! ; 1/EPS1 - 1
   REAL(wp), PUBLIC            ::   grav  = 9.80665_wp     !: gravity                            [m/s2]   
   
   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early
                                            !  output can be collected; do not change
                                            
   INTEGER       ::   nn_date0_hour=6  !: initial hour of the calendar da 
   
   PUBLIC   fld_read_rpn   ! called by sbc... modules

CONTAINS

   SUBROUTINE fld_read_rpn( kt, kn_fsbc, sd, GEM_rpn_list )
      implicit none
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_read_rpn  ***
      !!                   
      !! ** Purpose :   provide at each time step atmospheric variables
      !!                needed by sbc routine (core)
      !!                (UU,VV,TT, etc.) 
      !!
      !! ** Method  :   READ each input fields in STD files using RPN libs
      !!      and intepolate it to the model time-step.
      !!      Forcing files must be listed in "liste_inputfiles",located in
      !!      execution directory
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )               ::   kt        ! ocean time step
      INTEGER  , INTENT(in   )               ::   kn_fsbc   ! sbc computation period (in time step) 
      
   
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd        ! input field related variables
     
      !!
      INTEGER  ::   jf         ! dummy indices
      INTEGER  ::   ios        ! namelist error control
      
      character*512 GEM_rpn_list
      !!---------------------------------------------------------------------

      IF( kt == 1 ) THEN
      
        DO jf = 1, SIZE( sd )                     !    LOOP OVER FIELD    !
          sd(jf)%fnow(:,:,:)   = 0.0
          sd(jf)%fdta(:,:,:,:) = 0.0 ! tricky but working
        ENDDO
        
        kt_sbc=kt-1	
        GEM_rpn_list = GEM_rpn_list//'.txt'
        call init_atm (sd, GEM_rpn_list) ! get the initial forcing for step 0
        
      ENDIF

      kt_sbc=kt
      call linear_tint (sd) ! get the forcing for current time

       return
      END SUBROUTINE fld_read_rpn

      
      
      SUBROUTINE init_atm (sd, GEM_rpn_list)
      implicit none
      
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,dum
      character*16 datev
      character*512 GEM_rpn_list
      real*8  dayfrac,one,sid,rsid
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

      return
      END SUBROUTINE init_atm
      
      
   
      SUBROUTINE atm_openf(GEM_rpn_list)
      implicit none
      
      integer maxnfile
      parameter ( maxnfile=1000 )
      character*512 filename(maxnfile),fn,pwd,GEM_rpn_list
      integer  fnom,fstouv,fstlnk
      external fnom,fstouv,fstlnk
      integer err,err1,err2,i,cnt,unf,wkoffit
      nfile_std = 0
      cnt       = 0
      unf       = 0
      END SUBROUTINE atm_openf     
         
      
      

      SUBROUTINE linear_tint (sd)
      implicit none
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd     
      INTEGER  ::   jf         ! dummy indices

      character*16 datev,daten
      integer yy,mo,dd,hh,mm,ss,dum,datm
      real*8  tx,dayfrac,one,sid,rsid
      real(wp) ::  b
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)
          
      return

      END SUBROUTINE linear_tint



 
      
      
      
      
      SUBROUTINE atm_getdata (datev,put,get,sd)
      implicit none
      character*(*) datev
      logical put,get
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,datm,dum,jf
      tforc_1        = tforc_2
      return
      END SUBROUTINE atm_getdata

      
      
      
      SUBROUTINE get_col_position (lat_col,lon_col,i_col,j_col)
      implicit none

      real (kind=dbl_kind), intent(in):: lat_col, lon_col
      INTEGER, INTENT(out) :: i_col, j_col     
      INTEGER, DIMENSION(2) :: ind_col 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: wrk_col   
      
      allocate(wrk_col(jpi,jpj))   
      
      wrk_col = ((lat_rpn - lat_col)**2d0 + (lon_rpn - lon_col)**2d0)**0.5d0
      ind_col = minloc(wrk_col)
      i_col = ind_col(1)
      j_col = ind_col(2)
      print *, 'POSITION: ', lon_col, lat_col, lat_rpn(i_col,j_col), lon_rpn(i_col,j_col) 
      return
      END SUBROUTINE get_col_position      
      


   !!==============================================================================

END MODULE gemdrv_read_rpn
