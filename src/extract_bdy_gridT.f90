program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf coordinate file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridT files over specified period
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
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch, &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V, &
&                   data_suffix_ssh, data_suffix_bsf, data_suffix_ice,          &
&                   file_data_mask, file_data_zgr, file_data_hgr
CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_ssh, data_suffix_ice,  &
&                                       data_suffix_bsf, file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, myb, mxbt, glamt_ID, gphit_ID, &
&                                       e1t_ID, e2t_ID, nbit_ID, nbjt_ID, nbrt_ID, mtime, dimID_x, dimID_y,   &
&                                       mlon, mlat, mdeptht, kday, kmonth, kyear, kbdy, nfmt,        &
&                                       kt, kz,     lon_ID, lat_ID, deptht_ID, vosaline_ID, fidM,    &
&                                       votemper_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidT,  &
&                                       dimID_time_counter, dimID_deptht, time_ID, dimID_time, fidS, &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO,      &
&                                       depth_ID, ai, aj, bi, bj, kfmt, mx_REG, my_REG, tmask_REG_ID,&
&                                       tmask_GLO_ID, mx_tmp, my_tmp, fidMSKIN, fidZGR, e3t_GLO_ID,  &
&                                       mz_REG, fidMSKREG, dimID_z, mx_GLO, my_GLO, mz_GLO
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridT, file_bdy_gridT,   &
&                                       file_in_coord_REG, file_in_gridS, command_str, file_in_mask_REG
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: tmask_GLO, tmask_REG
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamt_bdy, gphit_bdy, e1t, e2t, nav_lon, nav_lat
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: votemper_GLO, vosaline_GLO, votemper_bdy, vosaline_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: deptht
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
REAL*8,ALLOCATABLE,DIMENSION(:,:,:)  :: e3t_GLO
LOGICAL                              :: existfile

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEt_ID, wgNWt_ID, wgSWt_ID, wgSEt_ID, angYt_ID, angXt_ID,  &
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
nn_bdy_eosmatch   =   1

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=bdy_data)
CLOSE(1)

!- bdy coordinates :
write(file_coord,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/coordinates_bdy_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coordinates_',a,'.nc')

write(file_in_mask_REG,104) TRIM(config_dir), TRIM(config)
104 FORMAT(a,'/mesh_mask_',a,'.nc')

eps=1.d-9

! name of interpolation weights file :
write(file_coeff,105) TRIM(config_dir), TRIM(config)
105 FORMAT(a,'/coeff_linear_',a,'.nc')

!=================================================================================
! 1a- Read GLOBAL mask :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_mask),0,fidMSKIN); call erreur(status,.TRUE.,"read mask input") 

status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z_GLO")
status = NF90_INQ_DIMID(fidMSKIN,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_GLO")
status = NF90_INQ_DIMID(fidMSKIN,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_GLO")

status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_z,len=mz_GLO); call erreur(status,.TRUE.,"inq_dim_z_GLO")
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_y,len=my_GLO); call erreur(status,.TRUE.,"inq_dim_y_GLO")
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_x,len=mx_GLO); call erreur(status,.TRUE.,"inq_dim_x_GLO")

ALLOCATE(  tmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask_GLO); call erreur(status,.TRUE.,"getvar_tmask_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1b- Read GLOBAL zgr (e3t) :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_zgr),0,fidZGR); call erreur(status,.TRUE.,"read e3t_GLO") 

ALLOCATE(  e3t_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidZGR,"e3t",e3t_GLO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGR,"e3t_0",e3t_GLO_ID)
call erreur(status,.TRUE.,"inq_e3t_GLO_ID")

status = NF90_GET_VAR(fidZGR,e3t_GLO_ID,e3t_GLO); call erreur(status,.TRUE.,"getvar_e3t_GLO")

status = NF90_CLOSE(fidZGR); call erreur(status,.TRUE.,"end read e3t_GLO")

!=================================================================================
! 1c- Read BDY coordinates
!=================================================================================

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbt",dimID_xbt) ; call erreur(status,.TRUE.,"inq_dimID_xbt")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbt,len=mxbt) ; call erreur(status,.TRUE.,"inq_dim_xbt")

ALLOCATE(  glamt_bdy(mxbt,myb)  ) 
ALLOCATE(  gphit_bdy(mxbt,myb)  ) 
ALLOCATE(  e1t(mxbt,myb)  ) 
ALLOCATE(  e2t(mxbt,myb)  ) 
ALLOCATE(  nbit(mxbt,myb)  ) 
ALLOCATE(  nbjt(mxbt,myb)  ) 
ALLOCATE(  nbrt(mxbt,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"e1t",e1t_ID)     ; call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidCOORD,"e2t",e2t_ID)     ; call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamt_bdy) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphit_bdy) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,e1t_ID,e1t)         ; call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidCOORD,e2t_ID,e2t)         ; call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidCOORD,nbit_ID,nbit)       ; call erreur(status,.TRUE.,"getvar_nbit")
status = NF90_GET_VAR(fidCOORD,nbjt_ID,nbjt)       ; call erreur(status,.TRUE.,"getvar_nbjt")
status = NF90_GET_VAR(fidCOORD,nbrt_ID,nbrt)       ; call erreur(status,.TRUE.,"getvar_nbrt")

status = NF90_CLOSE(fidCOORD) ; call erreur(status,.TRUE.,"close coordinate file")     

!=================================================================================
! 2- Read REGIONAL mask :
!=================================================================================

status = NF90_OPEN(TRIM(file_in_mask_REG),0,fidMSKREG); call erreur(status,.TRUE.,"read regional mask") 

status = NF90_INQ_DIMID(fidMSKREG,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z_REG")
status = NF90_INQ_DIMID(fidMSKREG,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_REG")
status = NF90_INQ_DIMID(fidMSKREG,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_REG")

status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_z,len=mz_REG); call erreur(status,.TRUE.,"inq_dim_z_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidMSKREG,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")

ALLOCATE(  tmask_REG(mx_REG,my_REG,mz_REG)  ) 

status = NF90_INQ_VARID(fidMSKREG,"tmask",tmask_REG_ID); call erreur(status,.TRUE.,"inq_tmask_REG_ID")

status = NF90_GET_VAR(fidMSKREG,tmask_REG_ID,tmask_REG); call erreur(status,.TRUE.,"getvar_tmask_REG")

status = NF90_CLOSE(fidMSKREG); call erreur(status,.TRUE.,"end read fidMSKREG")

!=====================================================================
! 3- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
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
! 4- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD_<data_suffix_T>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_<data_suffix_T>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <data_dir>/YYYY/<data_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <data_dir>/<data_prefix>_YYYY_MM_DD_<data_suffix_T>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <data_dir>/<data_prefix>_YYYY_MM_<data_suffix_T>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <data_dir>/<data_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <data_dir>/<data_prefix>_YYYY_MM.nc
291 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')          ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD_<data_suffix_T>.nc
292 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'_',a,'.nc')               ! <data_dir>/YYYY/<data_prefix>_YYYYMM_<data_suffix_T>.nc
293 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD.nc
294 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'.nc')                     ! <data_dir>/YYYY/<data_prefix>_YYYYMM.nc
295 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')                   ! <data_dir>/<data_prefix>_YYYYMMDD_<data_suffix_T>.nc
296 FORMAT(a,'/',a,'_',i4.4,i2.2,'_',a,'.nc')                        ! <data_dir>/<data_prefix>_YYYYMM_<data_suffix_T>.nc
297 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                         ! <data_dir>/<data_prefix>_YYYYMMDD.nc
298 FORMAT(a,'/',a,'_',i4.4,i2.2,'.nc')                              ! <data_dir>/<data_prefix>_YYYYMM.nc

kyear=nn_yeari
kmonth=1
DO kday=1,31

  ALLOCATE(list_fmt(16))
  list_fmt=(/191,192,193,194,195,196,197,198,291,292,293,294,295,296,297,298/)

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridT,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
        CASE(192)
          write(file_in_gridT,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
        CASE(193)
          write(file_in_gridT,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_gridT,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_gridT,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
        CASE(196) 
          write(file_in_gridT,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
        CASE(197) 
          write(file_in_gridT,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_gridT,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_gridT,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
        CASE(292)
          write(file_in_gridT,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
        CASE(293)
          write(file_in_gridT,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_gridT,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_gridT,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
        CASE(296) 
          write(file_in_gridT,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
        CASE(297) 
          write(file_in_gridT,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_gridT,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     inquire(file=file_in_gridT, exist=existfile)
     if ( existfile ) exit
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_gridT)
    status = NF90_OPEN(TRIM(file_in_gridT),0,fidT)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidT,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidT,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidT,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidT,"z",dimID_deptht)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidT,"depth",dimID_deptht)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidT,"deptht",dimID_deptht)
    call erreur(status,.TRUE.,"inq_dimID_deptht")

    status = NF90_INQUIRE_DIMENSION(fidT,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidT,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidT,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidT,dimID_deptht,len=mdeptht) ; call erreur(status,.TRUE.,"inq_dim_deptht")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( deptht(mdeptht) )

    status = NF90_INQ_VARID(fidT,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidT,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidT,"deptht",deptht_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"depth",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"nav_lev",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"z",depth_ID)
    call erreur(status,.TRUE.,"inq_deptht_ID")
        
    status = NF90_GET_VAR(fidT,lon_ID,nav_lon)                    ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidT,lat_ID,nav_lat)                    ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidT,deptht_ID,deptht)                  ; call erreur(status,.TRUE.,"getvar_deptht")

    status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")

    exit

  ELSEIF ( kday .eq. 31 ) THEN

    write(*,*) 'No temperature file found for first month of ', kyear
    write(*,*) '        >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!--

write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/BDY')
CALL system(TRIM(command_str))

!=================================================================================
! 3- Process all gridT files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridT,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(192)
          write(file_in_gridT,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
          write(file_in_gridS,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S) 
        CASE(193)
          write(file_in_gridT,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_gridT,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_gridT,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(196) 
          write(file_in_gridT,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
          write(file_in_gridS,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S)
        CASE(197) 
          write(file_in_gridT,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_gridT,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_gridT,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(292)
          write(file_in_gridT,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T) 
          write(file_in_gridS,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S) 
        CASE(293)
          write(file_in_gridT,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_gridT,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_gridT,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_T)
          write(file_in_gridS,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_S)
        CASE(296) 
          write(file_in_gridT,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_T)
          write(file_in_gridS,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_S)
        CASE(297) 
          write(file_in_gridT,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridS,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_gridT,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
          write(file_in_gridS,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridT, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          401 FORMAT(a,'/BDY/bdyT_tra_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridT,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          402 FORMAT(a,'/BDY/bdyT_tra_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridT,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridT  >>>> stop'
          stop
        endif

        ALLOCATE( votemper_GLO(mlon,mlat,mdeptht,mtime)  )
        ALLOCATE( vosaline_GLO(mlon,mlat,mdeptht,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input temperature :

        write(*,*) 'Reading temperature in ', TRIM(file_in_gridT)
        
        status = NF90_OPEN(TRIM(file_in_gridT),0,fidT)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidT,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidT,"votemper",votemper_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"toce",votemper_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"thetao",votemper_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"temperature",votemper_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"TEMPERATURE",votemper_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"TEMP",votemper_ID)
        call erreur(status,.TRUE.,"inq_votemper_ID")
        
        status = NF90_GET_VAR(fidT,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidT,votemper_ID,votemper_GLO)          ; call erreur(status,.TRUE.,"getvar_votemper")

        status = NF90_GET_ATT(fidT,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidT,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Read input salinity :

        write(*,*) 'Reading salinity in ', TRIM(file_in_gridS)
        
        status = NF90_OPEN(TRIM(file_in_gridS),0,fidS)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidS,"vosaline",vosaline_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidS,"soce",vosaline_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidS,"salinity",vosaline_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidS,"SALINITY",vosaline_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidS,"SALT",vosaline_ID)
        call erreur(status,.TRUE.,"inq_vosaline_ID")
        
        status = NF90_GET_VAR(fidS,vosaline_ID,vosaline_GLO)          ; call erreur(status,.TRUE.,"getvar_vosaline")

        status = NF90_CLOSE(fidS)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdeptht
        do l=1,mtime
          if ( .not. vosaline_GLO(i,j,k,l) .gt. 0.01 .or. .not. vosaline_GLO(i,j,k,l) .lt. 100.0 ) then
            vosaline_GLO(i,j,k,l) = 0.0
            votemper_GLO(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyT :
      
        ALLOCATE( votemper_bdy(mxbt,1,mdeptht,mtime)  )
        ALLOCATE( vosaline_bdy(mxbt,1,mdeptht,mtime)  )

        do kbdy=1,mxbt
        do kz=1,mdeptht

           aNW = wgNWt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz ) &
                                                  * e3t_GLO  ( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz )
           aSW = wgSWt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz ) &
                                                  * e3t_GLO  ( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz )
           aNE = wgNEt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz ) &
                                                  * e3t_GLO  ( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz )
           aSE = wgSEt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz ) &
                                                  * e3t_GLO  ( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz )
             
           aa = aNW + aSW + aNE + aSE

           if ( aa .gt. eps .and. zjSt(nbit(kbdy,1),nbjt(kbdy,1)) .gt. 1 ) then

             do kt=1,mtime
                votemper_bdy(kbdy,1,kz,kt) = (   votemper_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aNW   &
                &                              + votemper_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aSW   &
                &                              + votemper_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aNE   &
                &                              + votemper_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aSE ) &
                &                            * tmask_REG(nbit(kbdy,1),nbjt(kbdy,1),kz) / aa
                vosaline_bdy(kbdy,1,kz,kt) = (   vosaline_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aNW   &
                &                              + vosaline_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aSW   &
                &                              + vosaline_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aNE   &
                &                              + vosaline_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kz, kt ) * aSE ) &
                &                            * tmask_REG(nbit(kbdy,1),nbjt(kbdy,1),kz) / aa
             enddo

           else
 
             votemper_bdy(kbdy,1,kz,:) = 0.e0
             vosaline_bdy(kbdy,1,kz,:) = 0.e0
  
           endif

        enddo  !- kz
        enddo  !- kbdy

        !------------------------------------------------
        ! Convert to conservative temperature if needed :

        if ( nn_bdy_eosmatch .eq. 0 ) then
          write(*,*) 'Converting from EOS80 to TEOS10 ...'
          do kbdy=1,mxbt
          do kt=1,mtime
          do kz=1,mdeptht
            if ( vosaline_bdy(kbdy,1,kz,kt) .gt. 0.01 .and. vosaline_bdy(kbdy,1,kz,kt) .lt. 100.0 ) then
              vosaline_bdy(kbdy,1,kz,kt) = gsw_sa_from_sp( DBLE(vosaline_bdy(kbdy,1,kz,kt)), DBLE(deptht(kz)), DBLE(glamt_bdy(kbdy,1)), DBLE(gphit_bdy(kbdy,1)) )
              votemper_bdy(kbdy,1,kz,kt) = gsw_ct_from_pt( DBLE(vosaline_bdy(kbdy,1,kz,kt)), DBLE(votemper_bdy(kbdy,1,kz,kt)) )
            else
              vosaline_bdy(kbdy,1,kz,kt) = 0.0
              votemper_bdy(kbdy,1,kz,kt) = 0.0
            endif
          enddo
          enddo
          enddo
        elseif ( nn_bdy_eosmatch .ne. 1 ) then
          write(*,*) '~!@#$%^* Error: nn_bdy_eosmatch should be 0 or 1 >>>>> stop !!'
          stop
        endif

        !--------------------------------------
        ! Write BDY netcdf file

        write(*,*) 'Creating ', TRIM(file_bdy_gridT)
        status = NF90_CREATE(TRIM(file_bdy_gridT),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"deptht",mdeptht,dimID_deptht)                    ; call erreur(status,.TRUE.,"def_dimID_deptht")
        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbT",mxbT,dimID_xbT)                             ; call erreur(status,.TRUE.,"def_dimID_xbT")

        status = NF90_DEF_VAR(fidM,"votemper",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_deptht,dimID_time_counter/),votemper_ID)
        call erreur(status,.TRUE.,"def_var_votemper_ID")
        status = NF90_DEF_VAR(fidM,"vosaline",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_deptht,dimID_time_counter/),vosaline_ID)
        call erreur(status,.TRUE.,"def_var_vosaline_ID")
        status = NF90_DEF_VAR(fidM,"deptht",NF90_FLOAT,(/dimID_deptht/),deptht_ID)
        call erreur(status,.TRUE.,"def_var_deptht_ID")
        status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_xbT,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbrt_ID)
        call erreur(status,.TRUE.,"def_var_nbrt_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbT,dimID_yb/),nbjt_ID)
        call erreur(status,.TRUE.,"def_var_nbjt_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbT,dimID_yb/),nbit_ID)
        call erreur(status,.TRUE.,"def_var_nbit_ID")
       
        if ( nn_bdy_eosmatch .eq. 0 ) then 
          status = NF90_PUT_ATT(fidM,votemper_ID,"long_name","Conservative Temperature") ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,votemper_ID,"units","degC")                         ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","Absolute Salinity")        ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"units","g/kg")                         ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        else
          status = NF90_PUT_ATT(fidM,votemper_ID,"long_name","Temperature")              ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,votemper_ID,"units","degC")                         ; call erreur(status,.TRUE.,"put_att_votemper_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"long_name","Salinity")                 ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
          status = NF90_PUT_ATT(fidM,vosaline_ID,"units","psu")                          ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
        endif
        status = NF90_PUT_ATT(fidM,deptht_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_deptht_ID")
        status = NF90_PUT_ATT(fidM,deptht_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_deptht_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_gridT.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,votemper_ID,votemper_bdy) ; call erreur(status,.TRUE.,"var_votemper_ID")
        status = NF90_PUT_VAR(fidM,vosaline_ID,vosaline_bdy) ; call erreur(status,.TRUE.,"var_vosaline_ID")
        status = NF90_PUT_VAR(fidM,deptht_ID,deptht)         ; call erreur(status,.TRUE.,"var_deptht_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_bdy)   ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_bdy)   ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrt_ID,nbrt)             ; call erreur(status,.TRUE.,"var_nbrt_ID")
        status = NF90_PUT_VAR(fidM,nbjt_ID,nbjt)             ; call erreur(status,.TRUE.,"var_nbjt_ID")
        status = NF90_PUT_VAR(fidM,nbit_ID,nbit)             ; call erreur(status,.TRUE.,"var_nbit_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        DEALLOCATE( votemper_GLO, vosaline_GLO, time )
        DEALLOCATE( votemper_bdy, vosaline_bdy )

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
