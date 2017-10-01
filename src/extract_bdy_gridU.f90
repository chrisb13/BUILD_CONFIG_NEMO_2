program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to derive baroclinic velocities (U) along the BDY
!
! 0- Initialiartions
! 1a- Read GLOBAL mask
! 1b- Read GLOBAL zgr (e3u)
! 1c- Read GLOBAL hgr
! 1d- Read BDY coordinates
! 2- Read REGIONAL mask
! 3- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
! 4- Read input file dimensions in first existing file for specified time window
! 5- Process all gridU files over specified period
! 5a- Read input velocity
! 5b- Read input velocity
! 5c- Remove possible NaNs
! 5d- Express velocities onto the other grid
! 5e- Fill values on bdyU
! 5f- Write BDY netcdf file for baroclinic velocity
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch,    &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V,    &
&                   data_suffix_ssh, data_suffix_bsf, data_suffix_ice, sep1, sep2, &
&                   file_data_mask, file_data_zgr, file_data_hgr
CHARACTER(LEN=50)                    :: config, sep1, sep2
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_U, data_suffix_V, &
&                                       data_suffix_T, data_suffix_S, data_suffix_ssh, data_suffix_ice,  &
&                                       data_suffix_bsf, file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbu, myb, mxbu, glamt_ID, gphit_ID, &
&                                       e1u_ID, e2u_ID, nbiu_ID, nbju_ID, nbru_ID, mtime, dimID_x, dimID_y,   &
&                                       mlon, mlat, mdepthu, kday, kmonth, kyear, kbdy, nfmt, fidN,  &
&                                       kt, kz,     lon_ID, lat_ID, depthu_ID, vomecrty_ID,          &
&                                       vozocrtx_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidU,  &
&                                       dimID_time_counter, dimID_depthu, time_ID, dimID_time, fidV, &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO,      &
&                                       depth_ID, ai, aj, bi, bj, kfmt, mx_REG, my_REG, umask_REG_ID,&
&                                       umask_GLO_ID, mx_tmp, my_tmp, fidMSKIN, fidZGR, e3u_GLO_ID,  &
&                                       mz_REG, fidMSKREG, dimID_z, mx_GLO, my_GLO, mz_GLO, fidHGR,  &
&                                       e3u_REG_ID, e1u_GLO_ID, e1v_GLO_ID, im1, ip1, jm1, jp1
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridU, file_bdy_gridU,                   &
&                                       file_in_coord_REG, file_in_gridV, command_str,               &
&                                       file_bdy_gridU3d, file_in_mask_REG
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: umask_GLO, umask_REG
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbiu, nbju, nbru
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamu_bdy, gphiu_bdy, e1u_GLO, e1v_GLO, nav_lon, nav_lat
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: vozocrtx_GLO, vomecrty_GLO, vomecrty_GLO_U, vozocrtx_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: depthu
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
REAL*8,ALLOCATABLE,DIMENSION(:,:,:)  :: e3u_GLO, e3u_REG
LOGICAL                              :: existfile

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEu_ID, wgNWu_ID, wgSWu_ID, wgSEu_ID, angYu_ID, angXu_ID,  &
&                                         zkUt_ID, ziWu_ID, ziEu_ID, zjNu_ID, zjSu_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWu, ziEu, zjNu, zjSu
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYu, angXu
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEu, wgNWu, wgSWu, wgSEu
REAL*8                                 :: eps, aa, aSW, aNW, aSE, aNE, tmpbtp, tmpdep, uu, vv

!=================================================================================
! 0- Initialiartions
!=================================================================================

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

ALLOCATE(  umask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidMSKIN,"umask",umask_GLO_ID); call erreur(status,.TRUE.,"inq_umask_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,umask_GLO_ID,umask_GLO); call erreur(status,.TRUE.,"getvar_umask_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1b- Read GLOBAL zgr (e3u) :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_zgr),0,fidZGR); call erreur(status,.TRUE.,"read e3u_GLO") 

ALLOCATE(  e3u_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidZGR,"e3u",e3u_GLO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGR,"e3u_0",e3u_GLO_ID)
call erreur(status,.TRUE.,"inq_e3u_GLO_ID")

status = NF90_GET_VAR(fidZGR,e3u_GLO_ID,e3u_GLO); call erreur(status,.TRUE.,"getvar_e3u_GLO")

status = NF90_CLOSE(fidZGR); call erreur(status,.TRUE.,"end read e3u_GLO")

!=================================================================================
! 1c- Read GLOBAL hgr :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_hgr),0,fidHGR); call erreur(status,.TRUE.,"read hgr GLO") 

ALLOCATE(  e1u_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e1v_GLO(mx_GLO,my_GLO)  ) 

status = NF90_INQ_VARID(fidHGR,"e1u",e1u_GLO_ID) ; call erreur(status,.TRUE.,"inq_e1u_GLO_ID")
status = NF90_INQ_VARID(fidHGR,"e1v",e1v_GLO_ID) ; call erreur(status,.TRUE.,"inq_e1v_GLO_ID")

status = NF90_GET_VAR(fidHGR,e1u_GLO_ID,e1u_GLO); call erreur(status,.TRUE.,"getvar_e1u_GLO")
status = NF90_GET_VAR(fidHGR,e1v_GLO_ID,e1v_GLO); call erreur(status,.TRUE.,"getvar_e1v_GLO")

status = NF90_CLOSE(fidHGR); call erreur(status,.TRUE.,"end read hgr GLO")

!=================================================================================
! 1d- Read BDY coordinates
!=================================================================================

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbu",dimID_xbu) ; call erreur(status,.TRUE.,"inq_dimID_xbu")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbu,len=mxbu) ; call erreur(status,.TRUE.,"inq_dim_xbu")

ALLOCATE(  glamu_bdy(mxbu,myb)  ) 
ALLOCATE(  gphiu_bdy(mxbu,myb)  ) 
ALLOCATE(  nbiu(mxbu,myb)  ) 
ALLOCATE(  nbju(mxbu,myb)  ) 
ALLOCATE(  nbru(mxbu,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamu",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiu",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiu",nbiu_ID)   ; call erreur(status,.TRUE.,"inq_nbiu_ID")
status = NF90_INQ_VARID(fidCOORD,"nbju",nbju_ID)   ; call erreur(status,.TRUE.,"inq_nbju_ID")
status = NF90_INQ_VARID(fidCOORD,"nbru",nbru_ID)   ; call erreur(status,.TRUE.,"inq_nbru_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamu_bdy) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphiu_bdy) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,nbiu_ID,nbiu)       ; call erreur(status,.TRUE.,"getvar_nbiu")
status = NF90_GET_VAR(fidCOORD,nbju_ID,nbju)       ; call erreur(status,.TRUE.,"getvar_nbju")
status = NF90_GET_VAR(fidCOORD,nbru_ID,nbru)       ; call erreur(status,.TRUE.,"getvar_nbru")

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

ALLOCATE(  umask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  e3u_REG(mx_REG,my_REG,mz_REG)  ) 

status = NF90_INQ_VARID(fidMSKREG,"umask",umask_REG_ID); call erreur(status,.TRUE.,"inq_umask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"e3u",e3u_REG_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidMSKREG,"e3u_0",e3u_REG_ID)
call erreur(status,.TRUE.,"inq_e3u_REG_ID")

status = NF90_GET_VAR(fidMSKREG,umask_REG_ID,umask_REG); call erreur(status,.TRUE.,"getvar_umask_REG")
status = NF90_GET_VAR(fidMSKREG,e3u_REG_ID,e3u_REG); call erreur(status,.TRUE.,"getvar_e3u_REG")

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

ALLOCATE(  wgNEu(mx_REG,my_REG), wgNWu(mx_REG,my_REG), wgSWu(mx_REG,my_REG), wgSEu(mx_REG,my_REG) )
ALLOCATE(  ziWu(mx_REG,my_REG), ziEu(mx_REG,my_REG), zjNu(mx_REG,my_REG), zjSu(mx_REG,my_REG)  ) 
ALLOCATE(  angXu(mx_REG,my_REG), angYu(mx_REG,my_REG) )           
 
status = NF90_INQ_VARID(fidcoeff,"wgNEu",wgNEu_ID) ; call erreur(status,.TRUE.,"inq_wgNEu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWu",wgNWu_ID) ; call erreur(status,.TRUE.,"inq_wgNWu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWu",wgSWu_ID) ; call erreur(status,.TRUE.,"inq_wgSWu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEu",wgSEu_ID) ; call erreur(status,.TRUE.,"inq_wgSEu_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNu",zjNu_ID)   ; call erreur(status,.TRUE.,"inq_zjNu_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSu",zjSu_ID)   ; call erreur(status,.TRUE.,"inq_zjSu_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEu",ziEu_ID)   ; call erreur(status,.TRUE.,"inq_ziEu_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWu",ziWu_ID)   ; call erreur(status,.TRUE.,"inq_ziWu_ID")
status = NF90_INQ_VARID(fidcoeff,"angXu",angXu_ID) ; call erreur(status,.TRUE.,"inq_angXu_ID")
status = NF90_INQ_VARID(fidcoeff,"angYu",angYu_ID) ; call erreur(status,.TRUE.,"inq_angYu_ID")
               
status = NF90_GET_VAR(fidcoeff,wgNEu_ID,wgNEu)     ; call erreur(status,.TRUE.,"getvar_wgNEu")
status = NF90_GET_VAR(fidcoeff,wgNWu_ID,wgNWu)     ; call erreur(status,.TRUE.,"getvar_wgNWu")
status = NF90_GET_VAR(fidcoeff,wgSWu_ID,wgSWu)     ; call erreur(status,.TRUE.,"getvar_wgSWu")
status = NF90_GET_VAR(fidcoeff,wgSEu_ID,wgSEu)     ; call erreur(status,.TRUE.,"getvar_wgSEu")
status = NF90_GET_VAR(fidcoeff,zjNu_ID,zjNu)       ; call erreur(status,.TRUE.,"getvar_zjNu")
status = NF90_GET_VAR(fidcoeff,zjSu_ID,zjSu)       ; call erreur(status,.TRUE.,"getvar_zjSu")
status = NF90_GET_VAR(fidcoeff,ziEu_ID,ziEu)       ; call erreur(status,.TRUE.,"getvar_ziEu")
status = NF90_GET_VAR(fidcoeff,ziWu_ID,ziWu)       ; call erreur(status,.TRUE.,"getvar_ziWu")
status = NF90_GET_VAR(fidcoeff,angXu_ID,angXu)     ; call erreur(status,.TRUE.,"getvar_angXu")
status = NF90_GET_VAR(fidcoeff,angXu_ID,angXu)     ; call erreur(status,.TRUE.,"getvar_angXu")
                            
status = NF90_CLOSE(fidcoeff) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')  ! <data_dir>/YYYY/<data_prefix>YYYY<sep1>MM<sep2>DD<data_suffix>.nc  
192 FORMAT(a,'/',i4.4,'/',a,i4.4,a,i2.2,a,'.nc')         ! <data_dir>/YYYY/<data_prefix>YYYY<sep1>MM<data_suffix>.nc  
193 FORMAT(a,'/',a,i4.4,a,i2.2,a,i2.2,a,'.nc')           ! <data_dir>/<data_prefix>YYYY<sep1>MM<sep2>DD<data_suffix>.nc  
194 FORMAT(a,'/',a,i4.4,a,i2.2,a,'.nc')                  ! <data_dir>/<data_prefix>YYYY<sep1>MM<data_suffix>.nc 

ALLOCATE(list_fmt(4))
list_fmt=(/191,192,193,194/)

kyear=nn_yeari
kmonth=1
DO kday=1,31

  do kfmt=1,size(list_fmt)
     nfmt=list_fmt(kfmt)
     SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
        CASE(194)
          write(file_in_gridU,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_gridU)
     inquire(file=file_in_gridU, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading U input dimensions in ', TRIM(file_in_gridU)
    status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidU,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidU,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidU,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidU,"z",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depth",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depthu",dimID_depthu)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"deptht",dimID_depthu)
    call erreur(status,.TRUE.,"inq_dimID_depthu")

    status = NF90_INQUIRE_DIMENSION(fidU,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_depthu,len=mdepthu) ; call erreur(status,.TRUE.,"inq_dim_depthu")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( depthu(mdepthu) )

    status = NF90_INQ_VARID(fidU,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidU,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidU,"depthu",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"depth",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"deptht",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"nav_lev",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"z",depth_ID)
    call erreur(status,.TRUE.,"inq_depthu_ID")
        
    status = NF90_GET_VAR(fidU,lon_ID,nav_lon)                   ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidU,lat_ID,nav_lat)                   ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidU,depth_ID,depthu)                  ; call erreur(status,.TRUE.,"getvar_depthu")

    status = NF90_CLOSE(fidU)                                    ; call erreur(status,.TRUE.,"fin_lecture")

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
! 5- Process all gridU files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
          write(file_in_gridV,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U) 
          write(file_in_gridV,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_U)
          write(file_in_gridV,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_V)
        CASE(194)
          write(file_in_gridU,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_U)
          write(file_in_gridV,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_V)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridU, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          501 FORMAT(a,'/BDY/bdyU_u3d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU3d,501) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          502 FORMAT(a,'/BDY/bdyU_u3d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU3d,502) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridU3d and file_bdy_gridU3d  >>>> stop'
          stop
        endif

        ALLOCATE( vozocrtx_GLO(mlon,mlat,mdepthu,mtime)  )
        ALLOCATE( vomecrty_GLO(mlon,mlat,mdepthu,mtime)  )
        ALLOCATE( vomecrty_GLO_U(mlon,mlat,mdepthu,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! 5a- Read input velocity :

        write(*,*) 'Reading U in ', TRIM(file_in_gridU)
        
        status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)                ; call erreur(status,.TRUE.,"read U") 
        
        status = NF90_INQ_VARID(fidU,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidU,"vozocrtx",vozocrtx_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"uoce",vozocrtx_ID)
        call erreur(status,.TRUE.,"inq_vozocrtx_ID")
        
        status = NF90_GET_VAR(fidU,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidU,vozocrtx_ID,vozocrtx_GLO)          ; call erreur(status,.TRUE.,"getvar_vozocrtx")

        status = NF90_GET_ATT(fidU,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidU,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidU)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! 5b- Read input velocity :

        write(*,*) 'Reading V in ', TRIM(file_in_gridV)
        
        status = NF90_OPEN(TRIM(file_in_gridV),0,fidV)                ; call erreur(status,.TRUE.,"read V") 
        
        status = NF90_INQ_VARID(fidV,"vomecrty",vomecrty_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidV,"voce",vomecrty_ID)
        call erreur(status,.TRUE.,"inq_vomecrty_ID")
        
        status = NF90_GET_VAR(fidV,vomecrty_ID,vomecrty_GLO)          ; call erreur(status,.TRUE.,"getvar_vomecrty")

        status = NF90_CLOSE(fidV)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! 5c- Remove possible NaNs :

        write(*,*) 'Remove NaNs'
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdepthu
        do l=1,mtime
          if ( .not. abs(vomecrty_GLO(i,j,k,l)) .lt. 100.0 ) then
            vomecrty_GLO(i,j,k,l) = 0.0
          endif
          if ( .not. abs(vozocrtx_GLO(i,j,k,l)) .lt. 100.0 ) then
            vozocrtx_GLO(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! 5d- Express velocities onto the other grid :

        write(*,*) 'Express V on gridU'

        do i=1,mlon
        do j=1,mlat

          ip1=MIN(i+1,mlon)
          im1=MAX(i-1,1)
          jp1=MIN(j+1,mlat)
          jm1=MAX(j-1,1)

          vomecrty_GLO_U(i,j,:,:) = (   vomecrty_GLO(i  ,j  ,:,:) * e1v_GLO(i  ,j  ) &
          &                           + vomecrty_GLO(ip1,j  ,:,:) * e1v_GLO(ip1,j  ) &
          &                           + vomecrty_GLO(i  ,jm1,:,:) * e1v_GLO(i  ,jm1) &
          &                           + vomecrty_GLO(ip1,jm1,:,:) * e1v_GLO(ip1,jm1) ) / (4*e1u_GLO(i,j))

        enddo
        enddo

        DEALLOCATE( vomecrty_GLO )

        !---------------------------------------
        ! 5e- Fill values on bdyU :
      
        write(*,*) 'Interpolating then rotating...'

        ALLOCATE( vozocrtx_bdy(mxbu,1,mdepthu,mtime)  )

        do kbdy=1,mxbu

          do kz=1,mdepthu

            aNW = wgNWu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz ) &
                                                   * e3u_GLO  ( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz )
            aSW = wgSWu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz ) &
                                                   * e3u_GLO  ( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz )
            aNE = wgNEu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz ) &
                                                   * e3u_GLO  ( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz )
            aSE = wgSEu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz ) &
                                                   * e3u_GLO  ( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz )
             
            aa = aNW + aSW + aNE + aSE

            if ( aa .gt. eps .and. zjSu(nbiu(kbdy,1),nbju(kbdy,1)) .gt. 1 ) then

              do kt=1,mtime

                uu = (   vozocrtx_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aNW   &
                &      + vozocrtx_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aSW   &
                &      + vozocrtx_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aNE   &
                &      + vozocrtx_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aSE ) &
                &    * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),kz) / aa

                vv = (   vomecrty_GLO_U( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aNW   &
                &      + vomecrty_GLO_U( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aSW   &
                &      + vomecrty_GLO_U( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aNE   &
                &      + vomecrty_GLO_U( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kz, kt ) * aSE ) &
                &    * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),kz) / aa

                vozocrtx_bdy(kbdy,1,kz,kt) = uu * cos(angXu(nbiu(kbdy,1),nbju(kbdy,1))) + vv * sin(angYu(nbiu(kbdy,1),nbju(kbdy,1)))

              enddo !- kt

            else
 
              vozocrtx_bdy(kbdy,1,kz,:) = 0.e0
  
            endif

          enddo !- kz

          !- Removing barotropic component :
          !  NB: here no need to account for vvl because ubcl = u3d - sum_k(u3d*e3u)/sum_k(e3u)
          !                                                   = u3d - sum_k(u3d*e3u_0)/sum_k(e3u_0)
          !                                      because e3u = e3u_0 * ( 1 + ssh_u / H )
          do kt=1,mtime
            tmpdep = 0.e0
            tmpbtp = 0.e0
            do kz=1,mdepthu                   
               tmpbtp = tmpbtp + vozocrtx_bdy(kbdy,1,kz,kt) * e3u_REG(nbiu(kbdy,1),nbju(kbdy,1),kz) * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),kz)
               tmpdep = tmpdep +                              e3u_REG(nbiu(kbdy,1),nbju(kbdy,1),kz) * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),kz)
            enddo
            if ( tmpdep .gt. eps ) tmpbtp = tmpbtp / tmpdep
            do kz=1,mdepthu
              vozocrtx_bdy(kbdy,1,kz,kt) = ( vozocrtx_bdy(kbdy,1,kz,kt) - tmpbtp ) * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),kz)
            enddo
          enddo

        enddo !- kbdy

        DEALLOCATE( vozocrtx_GLO, vomecrty_GLO_U )

        !--------------------------------------
        ! 5f- Write BDY netcdf file for baroclinic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridU3d)
        status = NF90_CREATE(TRIM(file_bdy_gridU3d),NF90_NOCLOBBER,fidN) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidN,"depthu",mdepthu,dimID_depthu)                    ; call erreur(status,.TRUE.,"def_dimID_depthu")
        status = NF90_DEF_DIM(fidN,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidN,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidN,"xbU",mxbU,dimID_xbU)                             ; call erreur(status,.TRUE.,"def_dimID_xbU")

        status = NF90_DEF_VAR(fidN,"vozocrtx",NF90_FLOAT,(/dimID_xbU,dimID_yb,dimID_depthu,dimID_time_counter/),vozocrtx_ID)
        call erreur(status,.TRUE.,"def_var_vozocrtx_ID")
        status = NF90_DEF_VAR(fidN,"depthu",NF90_FLOAT,(/dimID_depthu/),depthu_ID)
        call erreur(status,.TRUE.,"def_var_depthu_ID")
        status = NF90_DEF_VAR(fidN,"nav_lat",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidN,"nav_lon",NF90_FLOAT,(/dimID_xbU,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidN,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidN,"nbrdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbru_ID)
        call erreur(status,.TRUE.,"def_var_nbru_ID")
        status = NF90_DEF_VAR(fidN,"nbjdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbju_ID)
        call erreur(status,.TRUE.,"def_var_nbju_ID")
        status = NF90_DEF_VAR(fidN,"nbidta",NF90_INT,(/dimID_xbU,dimID_yb/),nbiu_ID)
        call erreur(status,.TRUE.,"def_var_nbiu_ID")
       
        status = NF90_PUT_ATT(fidN,vozocrtx_ID,"long_name","Baroclinic Velocity")   ; call erreur(status,.TRUE.,"put_att_vozocrtx_ID")
        status = NF90_PUT_ATT(fidN,vozocrtx_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vozocrtx_ID")
        status = NF90_PUT_ATT(fidN,depthu_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_depthu_ID")
        status = NF90_PUT_ATT(fidN,depthu_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_depthu_ID")
        status = NF90_PUT_ATT(fidN,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidN,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,nbru_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidN,nbru_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidN,nbju_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidN,nbju_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidN,nbiu_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        status = NF90_PUT_ATT(fidN,nbiu_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"history","Created using extract_bdy_gridU.f90")
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidN) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidN,vozocrtx_ID,vozocrtx_bdy) ; call erreur(status,.TRUE.,"var_vozocrtx_ID")
        status = NF90_PUT_VAR(fidN,depthu_ID,depthu)         ; call erreur(status,.TRUE.,"var_depthu_ID")
        status = NF90_PUT_VAR(fidN,nav_lat_ID,gphiu_bdy)     ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidN,nav_lon_ID,glamu_bdy)     ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidN,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidN,nbru_ID,nbru)             ; call erreur(status,.TRUE.,"var_nbru_ID")
        status = NF90_PUT_VAR(fidN,nbju_ID,nbju)             ; call erreur(status,.TRUE.,"var_nbju_ID")
        status = NF90_PUT_VAR(fidN,nbiu_ID,nbiu)             ; call erreur(status,.TRUE.,"var_nbiu_ID")
        
        status = NF90_CLOSE(fidN) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        DEALLOCATE( vozocrtx_bdy, time )

        !--
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          write(*,*) 'Looking for next existing day in this month/year'
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
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
