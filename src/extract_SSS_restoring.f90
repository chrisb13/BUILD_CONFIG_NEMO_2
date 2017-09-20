program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf files with SSS (used for restoring)
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read REGIONAL mask
! 3- Read input file dimensions in first existing file for specified time window
!    (can be 2d or 3d salinity)
! 4- Process all gridT files over specified period
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /sss_resto/ nn_sss_yeari, nn_sss_yearf, sss_dir, sss_prefix, nn_sss_eosmatch, &
&                   sss_suffix, file_sss_mask
CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, sss_dir, sss_prefix, sss_suffix, file_sss_mask
INTEGER                              :: nn_sss_yeari, nn_sss_yearf, nn_sss_eosmatch

INTEGER                              :: status, mtime, dimID_x, dimID_y, mlon, mlat, kday, kmonth, kyear,   &
&                                       vosaline_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidSSS, mz,   &
&                                       dimID_time_counter, time_ID, dimID_time, fidS, kfmt, mx_tmp,        &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO, dimID_z,    &
&                                       nfmt, lon_ID, lat_ID, iREG, jREG, tmask_GLO_ID, ntest, ntest2,      &
&                                       fidMSKREG, mx_REG, my_REG, mz_REG, fidMSKIN, tmask_REG_ID, iii, jjj,&
&                                       kiter, rs, dij, im1, ip1, jm1, jp1, fidM, my_tmp
CHARACTER(LEN=100)                           :: calendar, time_units
CHARACTER(LEN=150)                           :: file_coord, file_in_SSS, file_REG_SSS, file_in_mask_REG, &
&                                               file_in_coord_REG, file_in_gridS, command_str
INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:)     :: list_fmt
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:)   :: tmask_GLO, tmask_REG, missing, tmp_missing
INTEGER(KIND=1),ALLOCATABLE,DIMENSION(:,:,:) :: tmask
REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:)      :: nav_lon, nav_lat, nav_lon_REG, nav_lat_REG
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)    :: vosaline_GLO, vosaline_REG, tmp_vosaline_REG
REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:,:,:)  :: vosaline_3d
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)        :: time
LOGICAL                                      :: existfile, ln_2d, iout

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEt_ID, wgNWt_ID, wgSWt_ID, wgSEt_ID, &
&                                         zkUt_ID, ziWt_ID, ziEt_ID, zjNt_ID, zjSt_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWt, ziEt, zjNt, zjSt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYt, angXt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEt, wgNWt, wgSWt, wgSEt
REAL*8                                 :: eps, aa, aSW, aNW, aSE, aNE

!=================================================================================
!- 0- Initialiartions
!=================================================================================

call gsw_saar_init (.true.)

write(*,*) 'Reading namelist parameters'

! Default values (replaced with namelist values if specified):
config_dir        = '.'
nn_sss_eosmatch   =   1

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=sss_resto)
CLOSE(1)

!- name of regional coordinates (input) :
write(file_in_coord_REG,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

! name of interpolation weights file :
write(file_coeff,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coeff_linear_',a,'.nc')

eps=1.d-9

!=================================================================================
! 1- Read REGIONAL mask :
!=================================================================================

write(*,*) 'Reading regional mask in ', TRIM(file_in_mask_REG)
status = NF90_OPEN(TRIM(file_in_mask_REG),0,fidMSKREG); call erreur(status,.TRUE.,"read regional mask") 
!-
status = NF90_INQ_DIMID(fidMSKREG,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z_REG")
status = NF90_INQ_DIMID(fidMSKREG,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_REG")
status = NF90_INQ_DIMID(fidMSKREG,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_REG")
!-
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_z,len=mz_REG); call erreur(status,.TRUE.,"inq_dim_z_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")
!-
ALLOCATE(  tmask(mx_REG,my_REG,mz_REG), tmask_REG(mx_REG,my_REG)  ) 
status = NF90_INQ_VARID(fidMSKREG,"tmask",tmask_REG_ID) ; call erreur(status,.TRUE.,"inq_tmask_REG_ID")
status = NF90_GET_VAR(fidMSKREG,tmask_REG_ID,tmask)     ; call erreur(status,.TRUE.,"getvar_tmask_REG")
tmask_REG(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
!-
status = NF90_CLOSE(fidMSKREG); call erreur(status,.TRUE.,"end read fidMSKREG")

!=====================================================================
! 2- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
!=====================================================================

write(*,*) 'Reading ', TRIM(file_coeff)

status = NF90_OPEN(TRIM(file_coeff),0,fidcoeff) ; call erreur(status,.TRUE.,"read weights") 
                                                  
status = NF90_INQ_DIMID(fidcoeff,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidcoeff,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_y,len=my_tmp); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_x,len=mx_tmp); call erreur(status,.TRUE.,"inq_dim_x_REG")
if ( mx_tmp .ne. mx_REG .or. my_tmp .ne. my_REG ) then
  write(*,*) '   ~!@#$%^* dimension mismatch >>>>>>> stop !!'
  stop
endif 

ALLOCATE(  wgNEt(mx_REG,my_REG), wgNWt(mx_REG,my_REG), wgSWt(mx_REG,my_REG), wgSEt(mx_REG,my_REG) )
ALLOCATE(  ziWt(mx_REG,my_REG), ziEt(mx_REG,my_REG), zjNt(mx_REG,my_REG), zjSt(mx_REG,my_REG)  ) 
            
status = NF90_INQ_VARID(fidcoeff,"wgNEt",wgNEt_ID) ; call erreur(status,.TRUE.,"inq_wgNEt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWt",wgNWt_ID) ; call erreur(status,.TRUE.,"inq_wgNWt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWt",wgSWt_ID) ; call erreur(status,.TRUE.,"inq_wgSWt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEt",wgSEt_ID) ; call erreur(status,.TRUE.,"inq_wgSEt_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNt",zjNt_ID)   ; call erreur(status,.TRUE.,"inq_zjNt_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSt",zjSt_ID)   ; call erreur(status,.TRUE.,"inq_zjSt_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEt",ziEt_ID)   ; call erreur(status,.TRUE.,"inq_ziEt_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWt",ziWt_ID)   ; call erreur(status,.TRUE.,"inq_ziWt_ID")
               
status = NF90_GET_VAR(fidcoeff,wgNEt_ID,wgNEt)     ; call erreur(status,.TRUE.,"getvar_wgNEt")
status = NF90_GET_VAR(fidcoeff,wgNWt_ID,wgNWt)     ; call erreur(status,.TRUE.,"getvar_wgNWt")
status = NF90_GET_VAR(fidcoeff,wgSWt_ID,wgSWt)     ; call erreur(status,.TRUE.,"getvar_wgSWt")
status = NF90_GET_VAR(fidcoeff,wgSEt_ID,wgSEt)     ; call erreur(status,.TRUE.,"getvar_wgSEt")
status = NF90_GET_VAR(fidcoeff,zjNt_ID,zjNt)       ; call erreur(status,.TRUE.,"getvar_zjNt")
status = NF90_GET_VAR(fidcoeff,zjSt_ID,zjSt)       ; call erreur(status,.TRUE.,"getvar_zjSt")
status = NF90_GET_VAR(fidcoeff,ziEt_ID,ziEt)       ; call erreur(status,.TRUE.,"getvar_ziEt")
status = NF90_GET_VAR(fidcoeff,ziWt_ID,ziWt)       ; call erreur(status,.TRUE.,"getvar_ziWt")
                            
status = NF90_CLOSE(fidcoeff) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 3- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <sss_dir>/YYYY/<sss_prefix>_YYYY_MM_DD_<sss_suffix>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <sss_dir>/YYYY/<sss_prefix>_YYYY_MM_<sss_suffix>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <sss_dir>/YYYY/<sss_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <sss_dir>/YYYY/<sss_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <sss_dir>/<sss_prefix>_YYYY_MM_DD_<sss_suffix>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <sss_dir>/<sss_prefix>_YYYY_MM_<sss_suffix>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <sss_dir>/<sss_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <sss_dir>/<sss_prefix>_YYYY_MM.nc
291 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')          ! <sss_dir>/YYYY/<sss_prefix>_YYYYMMDD_<sss_suffix>.nc
292 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'_',a,'.nc')               ! <sss_dir>/YYYY/<sss_prefix>_YYYYMM_<sss_suffix>.nc
293 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                ! <sss_dir>/YYYY/<sss_prefix>_YYYYMMDD.nc
294 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'.nc')                     ! <sss_dir>/YYYY/<sss_prefix>_YYYYMM.nc
295 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')                   ! <sss_dir>/<sss_prefix>_YYYYMMDD_<sss_suffix>.nc
296 FORMAT(a,'/',a,'_',i4.4,i2.2,'_',a,'.nc')                        ! <sss_dir>/<sss_prefix>_YYYYMM_<sss_suffix>.nc
297 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                         ! <sss_dir>/<sss_prefix>_YYYYMMDD.nc
298 FORMAT(a,'/',a,'_',i4.4,i2.2,'.nc')                              ! <sss_dir>/<sss_prefix>_YYYYMM.nc

kyear=nn_sss_yeari
kmonth=1
DO kday=1,31

  ALLOCATE(list_fmt(16))
  list_fmt=(/191,192,193,194,195,196,197,198,291,292,293,294,295,296,297,298/)

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_SSS,191) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(192)
          write(file_in_SSS,192) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix) 
        CASE(193)
          write(file_in_SSS,193) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_SSS,194) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_SSS,195) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(196) 
          write(file_in_SSS,196) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix)
        CASE(197) 
          write(file_in_SSS,197) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_SSS,198) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth
        CASE(291)
          write(file_in_SSS,291) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(292)
          write(file_in_SSS,292) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix) 
        CASE(293)
          write(file_in_SSS,293) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_SSS,294) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_SSS,295) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(296) 
          write(file_in_SSS,296) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix)
        CASE(297) 
          write(file_in_SSS,297) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_SSS,298) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     inquire(file=file_in_SSS, exist=existfile)
     if ( existfile ) exit
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_SSS)
    status = NF90_OPEN(TRIM(file_in_SSS),0,fidSSS)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidSSS,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidSSS,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidSSS,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    !- Figure out whether 2D or 3D salinity file :
    ln_2d = .false.
    status = NF90_INQ_DIMID(fidSSS,"z",dimID_z)
    if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidSSS,"depth",dimID_z)
    if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidSSS,"deptht",dimID_z)
    if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidSSS,"nav_lev",dimID_z)
    if ( status .ne. 0 ) ln_2d = .true.

    status = NF90_INQUIRE_DIMENSION(fidSSS,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidSSS,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidSSS,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    if ( .not. ln_2d ) then
      status = NF90_INQUIRE_DIMENSION(fidSSS,dimID_z,len=mz)
      call erreur(status,.TRUE.,"inq_dim_z")
      if ( mz .eq. 1 ) then
        ln_2d = .true.
      else
        write(*,*) 'Reading sea surface salinity from 3-dimensional files'
      endif
    endif

    write(*,*) 'dimensions :', mlon, mlat

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )

    status = NF90_INQ_VARID(fidSSS,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSSS,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSSS,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidSSS,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSSS,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidSSS,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
        
    status = NF90_GET_VAR(fidSSS,lon_ID,nav_lon)                    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidSSS,lat_ID,nav_lat)                    ; call erreur(status,.TRUE.,"getvar_lat")

    status = NF90_CLOSE(fidSSS)                                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No SSS file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!- Read tmask in large-scale/global file:
status = NF90_OPEN(TRIM(file_sss_mask),0,fidMSKIN);     call erreur(status,.TRUE.,"read mask_GLO") 
status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"depth",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"deptht",dimID_z)
if ( status .ne. 0 ) status = NF90_INQ_DIMID(fidMSKIN,"nav_lev",dimID_z)
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_z,len=mz)
ALLOCATE(  tmask(mlon,mlat,mz), tmask_GLO(mlon,mlat)  ) 
status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")
status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask);     call erreur(status,.TRUE.,"getvar_tmask_GLO")
tmask_GLO(:,:)=tmask(:,:,1)
DEALLOCATE(tmask)
status = NF90_CLOSE(fidMSKIN);                          call erreur(status,.TRUE.,"end read mask_GLO")

!- create SSS directory :
write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/SSS')
CALL system(TRIM(command_str))

!=================================================================================
! 4- Process all gridT files over specified period
!=================================================================================

DO kyear=nn_sss_yeari,nn_sss_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_SSS,191) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(192)
          write(file_in_SSS,192) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix) 
        CASE(193)
          write(file_in_SSS,193) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_SSS,194) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_SSS,195) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(196) 
          write(file_in_SSS,196) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix)
        CASE(197) 
          write(file_in_SSS,197) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_SSS,198) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth
        CASE(291)
          write(file_in_SSS,291) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(292)
          write(file_in_SSS,292) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix) 
        CASE(293)
          write(file_in_SSS,293) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_SSS,294) TRIM(sss_dir), kyear, TRIM(sss_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_SSS,295) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday, TRIM(sss_suffix)
        CASE(296) 
          write(file_in_SSS,296) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, TRIM(sss_suffix)
        CASE(297) 
          write(file_in_SSS,297) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_SSS,298) TRIM(sss_dir), TRIM(sss_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_SSS, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          401 FORMAT(a,'/SSS/sss_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_REG_SSS,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          402 FORMAT(a,'/SSS/sss_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_REG_SSS,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_REG_SSS  >>>> stop'
          stop
        endif

        ALLOCATE( vosaline_GLO(mlon,mlat,mtime)  )
        if ( .not. ln_2d ) ALLOCATE( vosaline_3d(mlon,mlat,mz,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input SSS :

        write(*,*) 'Reading SSS in ', TRIM(file_in_SSS)
        
        status = NF90_OPEN(TRIM(file_in_SSS),0,fidSSS)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidSSS,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidSSS,"vosaline",vosaline_ID)          ; call erreur(status,.TRUE.,"inq_vosaline_ID")
        
        status = NF90_GET_VAR(fidSSS,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        if ( ln_2d ) then
          status = NF90_GET_VAR(fidSSS,vosaline_ID,vosaline_GLO)        ; call erreur(status,.TRUE.,"getvar_vosaline")
        else
          status = NF90_GET_VAR(fidSSS,vosaline_ID,vosaline_3d)         ; call erreur(status,.TRUE.,"getvar_vosaline")
          vosaline_GLO(:,:,:) = vosaline_3d(:,:,1,:)
          DEALLOCATE( vosaline_3d )
        endif

        status = NF90_GET_ATT(fidSSS,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidSSS,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidSSS)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do l=1,mtime
          if ( .not. ( vosaline_GLO(i,j,l) .lt. 1.0d2 .and. vosaline_GLO(i,j,l) .gt. 1.0d-1 ) ) then
            vosaline_GLO(i,j,l) = 0.d0
          endif
        enddo
        enddo
        enddo

        !-------------------------------------------------
        ! convert to conservative temperature if needed :
        if ( nn_sss_eosmatch .eq. 0 ) then
          write(*,*) 'Converting from EOS80 to TEOS10 ...'
          do i=1,mlon
          do j=1,mlat
            do l=1,mtime
              if ( vosaline_GLO(i,j,l) .gt. 1.0d-1 .and. vosaline_GLO(i,j,l) .lt. 1.0d2 ) then
                vosaline_GLO(i,j,l) = gsw_sa_from_sp( DBLE(vosaline_GLO(i,j,l)), 1.d0, DBLE(nav_lon(i,j)), DBLE(nav_lat(i,j)) )
              else
                vosaline_GLO(i,j,l) = 0.d0
                tmask_GLO(i,j) = 0  ! force mask to 0 if meaningless value for at least 1 time step
              endif
            enddo
          enddo
          enddo
        elseif ( nn_sss_eosmatch .ne. 1 ) then
          write(*,*) '~!@#$%^* Error: nn_sss_eosmatch should be 0 or 1 >>>>> stop !!'
          stop
        endif
        
        !----------------------------------------------------------------------
        ! Just extract where ocean points on both global and regional grids :
        write(*,*) 'start extraction...'
        ALLOCATE( missing(mx_REG,my_REG) )
        ALLOCATE( tmp_missing(mx_REG,my_REG) )
        ALLOCATE( vosaline_REG(mx_REG,my_REG,mtime) )
        ALLOCATE( tmp_vosaline_REG(mx_REG,my_REG,mtime) )
        missing(:,:)=0
        vosaline_REG(:,:,:)=0.d0
        do iREG=1,mx_REG
        do jREG=1,my_REG

           aNW = tmask_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)
           aSW = tmask_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)
           aNE = tmask_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)
           aSE = tmask_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG)
           aa = aNW+aSW+aNE+aSE
           if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
             do l=1,mtime
               vosaline_REG(iREG,jREG,l)  =  (   vosaline_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG), l ) * aNW   &
               &                               + vosaline_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG), l ) * aSW   &
               &                               + vosaline_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG), l ) * aNE   &
               &                               + vosaline_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG), l ) * aSE ) &
               &                             * tmask_REG(iREG,jREG) / aa
             enddo
           elseif ( tmask_REG(iREG,jREG) .eq. 1 ) then
             missing(iREG,jREG) = 1
           endif

        enddo
        enddo

        ! Look for closest neighbours where we have missing values while tmask_REG=1:
        ntest2= NINT(sum(sum(1.0-FLOAT(tmask_REG(:,:)),2),1))
        do kiter=1,MAX(mx_REG/3,my_REG/3)
          ntest = NINT(sum(sum(FLOAT(missing),2),1))
          write(*,*) '  kiter = ', kiter
          write(*,*) '     nb of pts with missing value: ', ntest, ' for nb of land pts: ', ntest2
          if ( ntest .eq. 0 .or. ntest .eq. ntest2 ) exit
          tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
          tmp_missing(:,:)=missing(:,:)
          do iREG=1,mx_REG
          do jREG=1,my_REG
            if ( missing(iREG,jREG) .eq. 1 ) then
              iout=.FALSE.
              do rs=1,3,1
                  iii=iREG               ; jjj=jREG               
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=jREG               
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=jREG               
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=iREG               ; jjj=MIN(jREG+rs,my_REG)
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=iREG               ; jjj=MAX(jREG-rs,1)     
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=MIN(jREG+rs,my_REG)
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=MAX(jREG-rs,1)     
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=MIN(jREG+rs,my_REG)
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=MAX(jREG-rs,1)     
                  if ( tmask_REG(iii,jjj) .eq. 1 .and. missing(iii,jjj) .eq. 0 ) then
                    iout=.TRUE.
                    exit
                  endif
                enddo !- rs
              if (iout) then
                tmp_missing(iREG,jREG) = 0
                tmp_vosaline_REG(iREG,jREG,:) = vosaline_REG(iii,jjj,:)
                !exit
              elseif ( kiter .eq. MAX(mx_REG/3,my_REG/3) ) then
                write(*,953) iREG, jREG
                953 FORMAT(' >>> ERROR for point (',2I5,')')
                stop
              endif
            endif !-if ( missing(iREG,jREG) .eq. 1 )
          enddo !- jREG
          enddo !- iREG
          missing(:,:)=tmp_missing(:,:)
          vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
        enddo !- kiter
     
        !!- Smoothing (i.e. bi-linear interpolation) :
        !dij=INT(ai*0.5)
        !write(*,*) 'Smoothing over plus and minus :', dij, ' pts to avoid strong gradients'
        !tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
        !do iREG=1,mx_REG
        !do jREG=1,my_REG
        !    im1=MAX(iREG-dij,1) ; ip1=MIN(iREG+dij,mx_REG) 
        !    jm1=MAX(jREG-dij,1) ; jp1=MIN(jREG+dij,my_REG)
        !    if ( tmask_REG(iREG,jREG) .eq. 1 ) then 
        !      do l=1,mtime
        !        tmp_vosaline_REG(iREG,jREG,l) =   SUM( SUM( vosaline_REG(im1:ip1,jm1:jp1,l) * tmask_REG(im1:ip1,jm1:jp1), 2), 1) &
        !        &                               / SUM( SUM(                            1.0  * tmask_REG(im1:ip1,jm1:jp1), 2), 1)
        !      enddo
        !    else
        !      tmp_vosaline_REG(iREG,jREG,:) = 0.d0
        !    endif
        !enddo
        !enddo
        !vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
      
        !- "Drowning", i.e. put closest value everywhere on the mask file to avoid issue if namdom is slightly changed :
        !  We just repeat the previous methodology, but for masked points
        write(*,*) 'Drowning, i.e. fill all masked points with closest neighbour'
        missing(:,:)=NINT(1-FLOAT(tmask_REG(:,:)))
        ! Look for closest neighbours where we have missing values:
        do kiter=1,MAX(mx_REG,my_REG)
          ntest = NINT(sum(sum(FLOAT(missing),2),1))
          write(*,*) '  kiter = ', kiter
          write(*,*) '     remaining nb of masked points to fill: ', ntest
          if ( ntest .eq. 0 ) exit
          tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
          tmp_missing(:,:)=missing(:,:)
          do iREG=1,mx_REG
          do jREG=1,my_REG
            if ( missing(iREG,jREG) .eq. 1 ) then
              iout=.FALSE.
              do rs=1,3,1
                  iii=iREG               ; jjj=jREG               
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=jREG               
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=jREG               
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=iREG               ; jjj=MIN(jREG+rs,my_REG)
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=iREG               ; jjj=MAX(jREG-rs,1)     
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=MIN(jREG+rs,my_REG)
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MIN(iREG+rs,mx_REG); jjj=MAX(jREG-rs,1)     
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=MIN(jREG+rs,my_REG)
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
                  iii=MAX(iREG-rs,1)     ; jjj=MAX(jREG-rs,1)     
                  if ( missing(iii,jjj) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
              enddo !- rs
              if (iout) then
                tmp_missing(iREG,jREG) = 0
                tmp_vosaline_REG(iREG,jREG,:) = vosaline_REG(iii,jjj,:)
                !exit
              elseif ( kiter .eq. MAX(mx_REG,my_REG) ) then
                tmp_missing(iREG,jREG) = 0
                tmp_vosaline_REG(iREG,jREG,:) = 35.0000000000
                !exit
              endif
            endif !-if ( missing(iREG,jREG) .eq. 1 )
          enddo !- jREG
          enddo !- iREG
          missing(:,:)=tmp_missing(:,:)
          vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
        enddo !- kiter
      
        !--  
        DEALLOCATE( tmp_vosaline_REG, tmp_missing, missing )

        !--------------------------------------
        ! Write SSS netcdf file

        write(*,*) 'Creating ', TRIM(file_REG_SSS)
        status = NF90_CREATE(TRIM(file_REG_SSS),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create SSS file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x)                       ; call erreur(status,.TRUE.,"def_dimID_x")
        status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y)                       ; call erreur(status,.TRUE.,"def_dimID_y")
        
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_DOUBLE,(/dimID_time/),time_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"vosaline",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time/),vosaline_ID)
        call erreur(status,.TRUE.,"def_var_vosaline_ID")
        
        status = NF90_PUT_ATT(fidM,vosaline_ID,"associate","time_counter, z, y, x") ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        status = NF90_PUT_ATT(fidM,vosaline_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        status = NF90_PUT_ATT(fidM,vosaline_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        !status = NF90_PUT_ATT(fidM,vosaline_ID,"units","psu")                       ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        !status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","practical salinity")    ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        status = NF90_PUT_ATT(fidM,vosaline_ID,"units","g/kg")                      ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","absolute salinity")     ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"units",TRIM(time_units))                ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"calendar",TRIM(calendar))               ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"title","Time")                          ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"long_name","Time axis")                 ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"standard_name","time")                  ; call erreur(status,.TRUE.,"put_att_time_ID")
        status = NF90_PUT_ATT(fidM,time_ID,"axis","T")                              ; call erreur(status,.TRUE.,"put_att_time_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_SSS_restoring.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
        
        status = NF90_PUT_VAR(fidM,time_ID,time)             ; call erreur(status,.TRUE.,"var_time_ID")
        status = NF90_PUT_VAR(fidM,vosaline_ID,vosaline_REG) ; call erreur(status,.TRUE.,"var_vosaline_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close new SSS file")

        !--       
        DEALLOCATE( vosaline_GLO, time )
        DEALLOCATE( vosaline_REG )

        !--
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          write(*,*) 'Looking for next existing day in this month/year'
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          write(*,*) 'Only one file per month => switching to next month'
          exit
        else
          write(*,*) 'Keep in mind to decide how to end the main loop according to new file format  >>>> stop'
          stop
        endif

      ENDIF

    ENDDO !-kday
  ENDDO !- kmonth
ENDDO !- kyear

end program modif

!=============================================

SUBROUTINE erreur(iret, lstop, chaine)
  ! pour les messages d'erreur
  USE netcdf
  INTEGER, INTENT(in)                     :: iret
  LOGICAL, INTENT(in)                     :: lstop
  CHARACTER(LEN=*), INTENT(in)            :: chaine
  !
  CHARACTER(LEN=80)                       :: message
  !
  IF ( iret .NE. 0 ) THEN
    WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
    WRITE(*,*) 'ERREUR: ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
