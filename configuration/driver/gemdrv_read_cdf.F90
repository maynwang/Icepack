MODULE gemdrv_read_cdf
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
    USE gemdrv_read_rpn
    USE netcdf

   IMPLICIT NONE
   PUBLIC   




   INTEGER  ::   jpz = 1  
   INTEGER  ::   jpt = 25 


CONTAINS

   SUBROUTINE fld_read_cdf( kt, kn_fsbc, sd, GEM_cdf_list )
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
      
      character*512 GEM_cdf_list
      !!---------------------------------------------------------------------

      IF( kt == 1 ) THEN
      
        DO jf = 1, SIZE( sd )                     !    LOOP OVER FIELD    !
          sd(jf)%fnow(:,:,:)   = 0.0
          sd(jf)%fdta(:,:,:,:) = 0.0 ! tricky but working
        ENDDO
        
        kt_sbc=kt-1	
        call init_atm_cdf (sd, GEM_cdf_list) ! get the initial forcing for step 0
        
      ENDIF

      kt_sbc=kt
      call linear_tint_cdf(GEM_cdf_list,sd) ! get the forcing for current time

       return
      END SUBROUTINE fld_read_cdf

      
      
      SUBROUTINE init_atm_cdf (sd, GEM_cdf_list)
      implicit none
      
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,dum
      character*16 datev
      character*512 GEM_cdf_list
      real*8  dayfrac,one,sid,rsid
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

     !getting the initial date (Mod_runstrt_S) :
        write(Mod_runstrt_S(1:8),'(i8.8)') idate0
        Mod_runstrt_S(9:9)='.'
        write(Mod_runstrt_S(10:11),'(i2.2)') nh_gem_offs + nn_date0_hour
        Mod_runstrt_S(12:16)='0000 '  !Already assumed in opa


      ! The following is to get the date string (datev, which is the current file name) right
      
      dayfrac = dble(kt_sbc)*max(dt,3600.0)*rsid ! current timestep into seconds from sim start       
      call incdatsd  (datev,Mod_runstrt_S,dayfrac) !getting the date of current time
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1) !from date string to yy,mm,etc
      call pdfjdate2 (tforc_2,yy,mo,dd,hh,mm,ss) ! get current time (sec) into tforc_2
      ! Now we retrieve the first forcing data
      print *, 'Going into (initial) atm_get_data, datev, dayfrac : ', datev, dayfrac, dt
      call atm_getdata_cdf(GEM_cdf_list,datev,.false.,.true.,sd) ! get the date for current time (datev)
      
      print *, 'end atm_getdata initial'

      return
      END SUBROUTINE init_atm_cdf
      
      
      

      SUBROUTINE linear_tint_cdf(GEM_cdf_list,sd)
      implicit none
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd     
      INTEGER  ::   jf         ! dummy indices

      character*16 datev,daten
      character*512 GEM_cdf_list
      integer yy,mo,dd,hh,mm,ss,dum,datm
      real*8  tx,dayfrac,one,sid,rsid
      real(wp) ::  b
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)
          
      dayfrac = dble(kt_sbc)*max(dt,3600.0)*rsid 
      call incdatsd  (datev,Mod_runstrt_S,dayfrac)
      ! print *, 'datev, Mod_runstrt_S, dayfrac :', datev, Mod_runstrt_S, dayfrac
      if (datev.gt.current_atmf) then
      print *, 'Changing the forcing data'
      
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)      
      !if (yy < 2011)  then 
 	  dt_gem_atm = 3600.
      !else 
      !     dt_gem_atm = 3600.*3.0
      !endif      
         dayfrac = dble(dt_gem_atm)*rsid 
         call incdatsd (daten,current_atmf,dayfrac)
           ! print *, 'Going into atm_get_data, datev : ', datev     
         call atm_getdata_cdf (GEM_cdf_list,daten,.true.,.true.,sd)
      endif
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)
      call pdfjdate2 (tx,yy,mo,dd,hh,mm,ss)
      
      if (ln_gem_intrp) then
        b=(tx-tforc_1)/(tforc_2-tforc_1)
! Using averaged fields with date stamp corresponding to the end of
! the average serie
      else
        b=1.
      endif

      print *, 'interpolating with b = ', b
      DO jf = 1, SIZE( sd )
        call inter_field3 (b,sd(jf)%fdta(:,:,1,1), &
                             sd(jf)%fdta(:,:,1,2), &
                             sd(jf)%fnow(:,:,1),jpi,jpj)
      ENDDO

 1001 format (/' ----- ABORT ----- No ATM data valid at: ',a)
      return

      END SUBROUTINE linear_tint_cdf

      
      
      SUBROUTINE atm_getdata_cdf (GEM_cdf_list,datev,put,get,sd)
      implicit none
      ! In this function, we step the data in the sd array, and load the new forcing into
      ! from the time2 position.
      ! i.e. we do data 2 -> data1, and load new data into data2.

      character*(*) datev
      character*512 GEM_cdf_list
      logical put,get
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer yy,mo,dd,hh,mm,ss,datm,dum,jf
      
      tforc_1        = tforc_2
      DO jf = 1, SIZE( sd )     
        sd(jf)%fdta(:,:,1,1) = sd(jf)%fdta(:,:,1,2)
      ENDDO

      ! Load the data from the right date
      call atm_read_cdf (GEM_cdf_list,datev,sd)

      current_atmf = datev
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)
      call pdfjdate2 (tforc_2,yy,mo,dd,hh,mm,ss)

      return
      END SUBROUTINE atm_getdata_cdf



      SUBROUTINE atm_read_cdf (GEM_cdf_list,datev,sd)
      implicit none
      ! In this function, we open the cdfs in the folder
      ! corresponding to the datev, and load the data into the sd array

      character*(*) datev
      character*16 datev2
      character*512 GEM_cdf_list
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd
      integer datm, dum,ktgem
      real KNAMS,kprec
      parameter (KNAMS=0.514791)
      integer i,j,iavg,navg,nivr,nivt,nivm
      integer, dimension(2) :: lev_nul=(/-1,-1/), lev_wrk
      integer yy,mo,dd,hh,mm,ss
      INTEGER :: ni, nj,t
    
      real*8  dayfrac,one,sid,rsid
      character*16 datew
      character*4 year
      character*10 date00
      character (char_len_long) filename
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)
      real(wp) theta,pi,Tv


      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev,1)      

      ktgem =hh + 1
      ! print *, 'Finding the time level: ', hh, dt_gem_atm, ktgem
      ! print *, 'starting to write the cdf into arrays, at time level ', ktgem

      !Get a string with format YYYMMDD00      
      write(date00(1:10),10) yy,mo,dd
 10   format(i4.2,i2.2,i2.2,'00')

      !Get a string with format YYYY
      write(year(1:4),11) yy
 11   format(i4.2)

      !Get the netcdf file for the u wind component
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_u10.nc'
      !print *, 'filename: ', filename
      call update_variable_dimensions(filename)
      call readatm_fromCDF ( sd(1)%fdta (:,:,1,2),jpi,jpj,ktgem,'u_wind',filename,datev,KNAMS,0.)
	  !print *, 'UUOR sample : ', sd(1)%fdta(50,50,1,2)

      !Get the netcdf file for the v wind component
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_v10.nc'
      call readatm_fromCDF ( sd(2)%fdta (:,:,1,2),jpi,jpj,ktgem,'v_wind',filename,datev,KNAMS,0.)
	  !print *, 'VUOR sample : ', sd(2)%fdta(50,50,1,2)
	
      !Get the netcdf file for the humidity
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_q2.nc'
      call readatm_fromCDF ( sd(3)%fdta (:,:,1,2),jpi,jpj,ktgem,'qair'  ,filename,datev,1.0,0.)
	  !print *, 'HU sample : ', sd(2)%fdta(50,50,1,2)

      !Get the netcdf file for the short-wave radiations
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_qsw.nc'
      call readatm_fromCDF ( sd(4)%fdta (:,:,1,2),jpi,jpj,ktgem,'solar'  ,filename,datev,1.0  ,0.)
	  !print *, 'FB sample : ', sd(4)%fdta(50,50,1,2)
	  
      !Get the netcdf file for the long-wave radiations
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_qlw.nc'
      call readatm_fromCDF ( sd(5)%fdta (:,:,1,2),jpi,jpj,ktgem,'therm_rad'  ,filename,datev,1.0  ,0.)
	  !print *, 'FI sample : ', sd(5)%fdta(50,50,1,2)
	  
      !Get the netcdf file for the air temperature
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_t2.nc'
      call readatm_fromCDF ( sd(6)%fdta (:,:,1,2),jpi,jpj,ktgem,'tair'  ,filename,datev,1.0,0.)
	  !print *, 'TT sample : ', sd(6)%fdta(50,50,1,2)

      !Get the netcdf file for the precipitations
      filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_precip.nc'
      call readatm_fromCDF ( sd(7)%fdta (:,:,1,2),jpi,jpj,ktgem,'precip'  ,filename,datev,kprec,0.)
	  ! print *, 'PR sample : ', sd(7)%fdta(50,50,1,2) 
          !print *, 'filename: ', filename
          
      !print *, 'zlev_gem :', zlev_gem
      !Getting the sea level pressure
      if (zlev_gem.lt.0.) then
          !lev_wrk(:) = nivt
          
          !Get the netcdf file for the sea level pressure (therm level height)
          filename = trim(GEM_cdf_list)//'/'//trim(year)//'/'//trim(date00)//'/'//trim(date00)//'_slp.nc'
          call readatm_fromCDF (sd(11)%fdta(:,:,1,2),jpi,jpj,ktgem,'atmpres' ,filename,datev,100.0,0.)
	      !print *, 'PX sample : ', sd(11)%fdta(50,50,1,2)
	      
          !Get the netcdf file for the sea level pressure (mometum level height)
          sd(12)%fdta(:,:,1,2)=sd(11)%fdta(:,:,1,2)
        
          !Get the netcdf file for the sea level pressure (sea level height)
          sd(13)%fdta(:,:,1,2)=sd(11)%fdta(:,:,1,2)
         !        if(ln_apr_dyn) &
         !        call readfstatm (sd(15)%fdta(:,:,1,2),jpi,jpj,'PN'  ,datev,lev_nul, 100.0,0.,nivr )
         !   ... Thickness of first atm. model layer
         !   ... (hydrostatic relation approximating Tv_bar=Tv(layer=1))
          do j=1,jpj
          do i=1,jpi
             Tv= sd(6)%fdta(i,j,1,2) * (1.0D0 + DELTA*sd(3)%fdta(i,j,1,2))
             sd(9 )%fdta(i,j,1,2) = -RGASD*Tv/grav*log(sd(11)%fdta(i,j,1,2)/sd(13)%fdta(i,j,1,2))
             sd(10)%fdta(i,j,1,2) = -RGASD*Tv/grav*log(sd(12)%fdta(i,j,1,2)/sd(13)%fdta(i,j,1,2))
          enddo
          enddo

      else
        sd(9 )%fdta(:,:,1,2) = zlev_gem
        sd(10)%fdta(:,:,1,2) = zlev_gem
        sd(11)%fdta(:,:,1,2) = 100000.
        sd(12)%fdta(:,:,1,2) = 100000.
        sd(13)%fdta(:,:,1,2) = 100000.
        sd(15)%fdta(:,:,1,2) = 100000.
      endif


      do j=1,jpj 
      do i=1,jpi 
        sd(8)%fdta(i,j,1,2)=0.
        if (sd(6)%fdta(i,j,1,2).le.273.16) sd(8)%fdta(i,j,1,2)=sd(7)%fdta(i,j,1,2)
                ! Now using upper level, 8 corresponds to snow
      enddo
      enddo

       ! Converts sd(6)%fdta (:,:,1,2) to potential temperature
      sd(6)%fdta(:,:,1,2) = sd(6)%fdta(:,:,1,2)*    &
     &                (sd(13)%fdta(:,:,1,2)/sd(11)%fdta(:,:,1,2))**CAPPA
     
     
      !extraction of the lat and lon arrays      
      if (datev .eq. Mod_runstrt_S) then
            call prsdate   (2001,01,10,00,00,00,dum,datev2,2)
      	  !print *, 'DATEV2: ', datev2
          allocate(lat_rpn(jpi,jpj))
          allocate(lon_rpn(jpi,jpj))        
          call readlatlon( lat_rpn(:,:),jpi,jpj,'nav_lat',filename)
	   !   print *, 'LAT sample : ', lat_rpn(500,500)   
          call readlatlon( lon_rpn(:,:),jpi,jpj,'nav_lon',filename)
	    !  print *, 'LON sample : ', lon_rpn(500,500)     
      endif  
	  print *, 'END get DATA!!!!'
	
      return
      
      END SUBROUTINE atm_read_cdf

  FUNCTION netcdf_check(status, ncid) RESULT(ret)
    INTEGER, INTENT(in) :: status
    INTEGER, INTENT(in) :: ncid
    INTEGER :: ret
    INTEGER :: close_status

    ret = NF90_NOERR

    IF (status /= NF90_NOERR) THEN
      PRINT *, 'NetCDF Error: ', trim(nf90_strerror(status))
      IF (ncid /= -1) THEN
        close_status = nf90_close(ncid)
        IF (close_status /= NF90_NOERR) THEN
          PRINT *, 'Error while closing file:', trim(nf90_strerror(close_status))
        ELSE
          PRINT *, 'File closed due to error.'
        END IF
      END IF
      CALL exit(1)
    END IF

    RETURN
  END FUNCTION netcdf_check   
      
    SUBROUTINE readatm_fromCDF(f, jpi, jpj, t, varname, filename, dat, factm, facta)
   
    IMPLICIT NONE 

    CHARACTER(*) varname, filename, dat
    integer jpi, jpj
    real factm, facta
    REAL(wp), DIMENSION(jpi,jpj) :: f
    real(wp), DIMENSION(jpi,jpj,jpz,jpt) :: wrk

    INTEGER :: i, j, z, t
    INTEGER :: ncid, niid, njid, nzid, ntid, varid    ! IDs for netcdf file, dimensions, layer variable
    INTEGER :: status ! Variable for netcdf subroutine status

    PRINT *, 'Getting into CDF files for ', dat, t, filename, varname

    status = netcdf_check(nf90_open(path = filename, mode = nf90_nowrite, ncid = ncid), -1)
    status = netcdf_check(nf90_inq_dimid(ncid, "x", niid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "y", njid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "z", nzid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "time_counter", ntid), ncid)

    status = netcdf_check(nf90_inquire_dimension(ncid, niid, len = jpi), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, njid, len = jpj), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, nzid, len = jpz), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, ntid, len = jpt), ncid)

    status = netcdf_check(nf90_inq_varid(ncid, varname, varid), ncid)
    status = netcdf_check(nf90_get_var(ncid, varid, wrk), ncid)
    status = netcdf_check(nf90_close(ncid), ncid)

    z = 1
    DO i = 1, jpi
      DO j = 1, jpj
        f(i, j) = wrk(i, j, z, t)
      END DO
    END DO

    f = f * factm + facta
  END SUBROUTINE readatm_fromCDF


  SUBROUTINE readlatlon(f, jpi, jpj, varname, filename)
      implicit none

      character*(*) varname,filename
      integer jpi,jpj
      REAL(wp) , DIMENSION(jpi,jpj) :: f
      real(wp) , DIMENSION(jpi,jpj) :: wrk

      integer i,j
      integer ncid, niid, njid, nzid, ntid, varid       ! IDs for netcdf file, dimensions, layer variable
      INTEGER :: status                             ! Variable for netcdf subroutine status
    
    status = netcdf_check(nf90_open(path = filename, mode = nf90_nowrite, ncid = ncid), -1)
    status = netcdf_check(nf90_inq_dimid(ncid, "x", niid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "y", njid), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, niid, len = jpi), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, njid, len = jpj), ncid)
    status = netcdf_check(nf90_inq_varid(ncid, varname, varid), ncid)
    status = netcdf_check(nf90_get_var(ncid, varid, wrk), ncid)
    status = netcdf_check(nf90_close(ncid), ncid)

    DO i = 1, jpi
      DO j = 1, jpj
        f(i, j) = wrk(i, j)
      END DO
    END DO
  END SUBROUTINE readlatlon


  SUBROUTINE update_variable_dimensions(filename)
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: filename
    INTEGER :: jpi, jpj, jpz, jpt
    INTEGER :: ncid, niid, njid, nzid, ntid
    INTEGER :: status

    status = netcdf_check(nf90_open(path = filename, mode = nf90_nowrite, ncid = ncid), -1)
    status = netcdf_check(nf90_inq_dimid(ncid, "x", niid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "y", njid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "z", nzid), ncid)
    status = netcdf_check(nf90_inq_dimid(ncid, "time_counter", ntid), ncid)

    status = netcdf_check(nf90_inquire_dimension(ncid, niid, len = jpi), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, njid, len = jpj), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, nzid, len = jpz), ncid)
    status = netcdf_check(nf90_inquire_dimension(ncid, ntid, len = jpt), ncid)
    status = netcdf_check(nf90_close(ncid), ncid)
  END SUBROUTINE update_variable_dimensions

END MODULE gemdrv_read_cdf
