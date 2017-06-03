program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf coordinate file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all gridU files over specified period
!
! history : - Feb. 2017: version with namelist (N. Jourdain)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /bdy_data/ nn_yeari, nn_yearf, data_dir, data_prefix, nn_bdy_eosmatch, &
&                   data_suffix_T, data_suffix_S, data_suffix_U, data_suffix_V, &
&                   data_suffix_ssh, data_suffix_bsf, data_suffix_ice,          &
&                   file_data_mask, file_data_zgr, file_data_hgr
CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_U, data_suffix_V, &
&                                       data_suffix_T, data_suffix_S, data_suffix_ssh, data_suffix_ice,  &
&                                       data_suffix_bsf, file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbv, myb, mxbv, glamt_ID, gphit_ID, &
&                                       e1u_ID, e2u_ID, nbiv_ID, nbjv_ID, nbrv_ID, mtime, dimID_x, dimID_y,   &
&                                       mlon, mlat, mdepthv, kday, kmonth, kyear, kbdy, nfmt, fidN,  &
&                                       kt, kz,     lon_ID, lat_ID, depthv_ID, vozocrtx_ID, fidM,    &
&                                       vomecrty_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidU,  &
&                                       dimID_time_counter, dimID_depthv, time_ID, dimID_time, fidV, &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO,      &
&                                       depth_ID, ai, aj, bi, bj, kfmt, mx_REG, my_REG, vmask_REG_ID,&
&                                       vmask_GLO_ID, mx_tmp, my_tmp, fidMSKIN, fidZGR, e3v_GLO_ID,  &
&                                       mz_REG, fidMSKREG, dimID_z, mx_GLO, my_GLO, mz_GLO, fidHGR,  &
&                                       e3v_REG_ID, e2u_GLO_ID, e2v_GLO_ID, im1, ip1, jm1, jp1
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_gridU, file_bdy_gridV, file_bdy_gridV2d, &
&                                       file_in_coord_REG, file_in_gridV, command_str,               &
&                                       file_bdy_gridV3d, file_in_mask_REG
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: vmask_GLO, vmask_REG
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbiv, nbjv, nbrv
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamv_bdy, gphiv_bdy, e2u_GLO, e2v_GLO, nav_lon, nav_lat
REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:):: vomecrty_GLO, vozocrtx_GLO, vozocrtx_GLO_V, vomecrty_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: depthv
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
REAL*8,ALLOCATABLE,DIMENSION(:,:,:)  :: e3v_GLO, e3v_REG
LOGICAL                              :: existfile

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEv_ID, wgNWv_ID, wgSWv_ID, wgSEv_ID, angYv_ID, angXv_ID,  &
&                                         zkUt_ID, ziWv_ID, ziEv_ID, zjNv_ID, zjSv_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWv, ziEv, zjNv, zjSv
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYv, angXv
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEv, wgNWv, wgSWv, wgSEv
REAL*8                                 :: eps, aa, aSW, aNW, aSE, aNE, uu, vv, tmpdep, tmpbtp

!=================================================================================
!- 0- Initialiartions
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

ALLOCATE(  vmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidMSKIN,"vmask",vmask_GLO_ID); call erreur(status,.TRUE.,"inq_vmask_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,vmask_GLO_ID,vmask_GLO); call erreur(status,.TRUE.,"getvar_vmask_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1b- Read GLOBAL zgr (e3t) :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_zgr),0,fidZGR); call erreur(status,.TRUE.,"read e3v_GLO") 

ALLOCATE(  e3v_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidZGR,"e3v",e3v_GLO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGR,"e3v_0",e3v_GLO_ID)
call erreur(status,.TRUE.,"inq_e3v_GLO_ID")

status = NF90_GET_VAR(fidZGR,e3v_GLO_ID,e3v_GLO); call erreur(status,.TRUE.,"getvar_e3v_GLO")

status = NF90_CLOSE(fidZGR); call erreur(status,.TRUE.,"end read e3v_GLO")

!=================================================================================
! 1c- Read GLOBAL hgr :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_data_hgr),0,fidHGR); call erreur(status,.TRUE.,"read hgr GLO") 

ALLOCATE(  e2u_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e2v_GLO(mx_GLO,my_GLO)  ) 

status = NF90_INQ_VARID(fidHGR,"e2u",e2u_GLO_ID) ; call erreur(status,.TRUE.,"inq_e2u_GLO_ID")
status = NF90_INQ_VARID(fidHGR,"e2v",e2v_GLO_ID) ; call erreur(status,.TRUE.,"inq_e2v_GLO_ID")

status = NF90_GET_VAR(fidHGR,e2u_GLO_ID,e2u_GLO); call erreur(status,.TRUE.,"getvar_e2u_GLO")
status = NF90_GET_VAR(fidHGR,e2v_GLO_ID,e2v_GLO); call erreur(status,.TRUE.,"getvar_e2v_GLO")

status = NF90_CLOSE(fidHGR); call erreur(status,.TRUE.,"end read hgr GLO")

!=================================================================================
! 1d- Read BDY coordinates
!=================================================================================

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbv",dimID_xbv) ; call erreur(status,.TRUE.,"inq_dimID_xbv")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbv,len=mxbv) ; call erreur(status,.TRUE.,"inq_dim_xbv")

ALLOCATE(  glamv_bdy(mxbv,myb)  ) 
ALLOCATE(  gphiv_bdy(mxbv,myb)  ) 
ALLOCATE(  nbiv(mxbv,myb)  ) 
ALLOCATE(  nbjv(mxbv,myb)  ) 
ALLOCATE(  nbrv(mxbv,myb)  ) 

status = NF90_INQ_VARID(fidCOORD,"glamv",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiv",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiv",nbiv_ID)   ; call erreur(status,.TRUE.,"inq_nbiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjv",nbjv_ID)   ; call erreur(status,.TRUE.,"inq_nbjv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrv",nbrv_ID)   ; call erreur(status,.TRUE.,"inq_nbrv_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamv_bdy) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphiv_bdy) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidCOORD,nbiv_ID,nbiv)       ; call erreur(status,.TRUE.,"getvar_nbiv")
status = NF90_GET_VAR(fidCOORD,nbjv_ID,nbjv)       ; call erreur(status,.TRUE.,"getvar_nbjv")
status = NF90_GET_VAR(fidCOORD,nbrv_ID,nbrv)       ; call erreur(status,.TRUE.,"getvar_nbrv")

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

ALLOCATE(  vmask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  e3v_REG(mx_REG,my_REG,mz_REG)  ) 

status = NF90_INQ_VARID(fidMSKREG,"vmask",vmask_REG_ID); call erreur(status,.TRUE.,"inq_vmask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"e3v",e3v_REG_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidMSKREG,"e3v_0",e3v_REG_ID)
call erreur(status,.TRUE.,"inq_e3v_REG_ID")

status = NF90_GET_VAR(fidMSKREG,vmask_REG_ID,vmask_REG); call erreur(status,.TRUE.,"getvar_vmask_REG")
status = NF90_GET_VAR(fidMSKREG,e3v_REG_ID,e3v_REG); call erreur(status,.TRUE.,"getvar_e3v_REG")

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

ALLOCATE(  wgNEv(mx_REG,my_REG), wgNWv(mx_REG,my_REG), wgSWv(mx_REG,my_REG), wgSEv(mx_REG,my_REG) )
ALLOCATE(  ziWv(mx_REG,my_REG), ziEv(mx_REG,my_REG), zjNv(mx_REG,my_REG), zjSv(mx_REG,my_REG)  ) 
ALLOCATE(  angXv(mx_REG,my_REG), angYv(mx_REG,my_REG) )           
 
status = NF90_INQ_VARID(fidcoeff,"wgNEv",wgNEv_ID) ; call erreur(status,.TRUE.,"inq_wgNEv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWv",wgNWv_ID) ; call erreur(status,.TRUE.,"inq_wgNWv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWv",wgSWv_ID) ; call erreur(status,.TRUE.,"inq_wgSWv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEv",wgSEv_ID) ; call erreur(status,.TRUE.,"inq_wgSEv_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNv",zjNv_ID)   ; call erreur(status,.TRUE.,"inq_zjNv_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSv",zjSv_ID)   ; call erreur(status,.TRUE.,"inq_zjSv_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEv",ziEv_ID)   ; call erreur(status,.TRUE.,"inq_ziEv_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWv",ziWv_ID)   ; call erreur(status,.TRUE.,"inq_ziWv_ID")
status = NF90_INQ_VARID(fidcoeff,"angXv",angXv_ID) ; call erreur(status,.TRUE.,"inq_angXv_ID")
status = NF90_INQ_VARID(fidcoeff,"angYv",angYv_ID) ; call erreur(status,.TRUE.,"inq_angYv_ID")
               
status = NF90_GET_VAR(fidcoeff,wgNEv_ID,wgNEv)     ; call erreur(status,.TRUE.,"getvar_wgNEv")
status = NF90_GET_VAR(fidcoeff,wgNWv_ID,wgNWv)     ; call erreur(status,.TRUE.,"getvar_wgNWv")
status = NF90_GET_VAR(fidcoeff,wgSWv_ID,wgSWv)     ; call erreur(status,.TRUE.,"getvar_wgSWv")
status = NF90_GET_VAR(fidcoeff,wgSEv_ID,wgSEv)     ; call erreur(status,.TRUE.,"getvar_wgSEv")
status = NF90_GET_VAR(fidcoeff,zjNv_ID,zjNv)       ; call erreur(status,.TRUE.,"getvar_zjNv")
status = NF90_GET_VAR(fidcoeff,zjSv_ID,zjSv)       ; call erreur(status,.TRUE.,"getvar_zjSv")
status = NF90_GET_VAR(fidcoeff,ziEv_ID,ziEv)       ; call erreur(status,.TRUE.,"getvar_ziEv")
status = NF90_GET_VAR(fidcoeff,ziWv_ID,ziWv)       ; call erreur(status,.TRUE.,"getvar_ziWv")
status = NF90_GET_VAR(fidcoeff,angXv_ID,angXv)     ; call erreur(status,.TRUE.,"getvar_angXv")
status = NF90_GET_VAR(fidcoeff,angYv_ID,angYv)     ; call erreur(status,.TRUE.,"getvar_angYv")
                            
status = NF90_CLOSE(fidcoeff) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD_<data_suffix_U>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_<data_suffix_U>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <data_dir>/YYYY/<data_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <data_dir>/<data_prefix>_YYYY_MM_DD_<data_suffix_U>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <data_dir>/<data_prefix>_YYYY_MM_<data_suffix_U>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <data_dir>/<data_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <data_dir>/<data_prefix>_YYYY_MM.nc
291 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')          ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD_<data_suffix_U>.nc
292 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'_',a,'.nc')               ! <data_dir>/YYYY/<data_prefix>_YYYYMM_<data_suffix_U>.nc
293 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD.nc
294 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'.nc')                     ! <data_dir>/YYYY/<data_prefix>_YYYYMM.nc
295 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')                   ! <data_dir>/<data_prefix>_YYYYMMDD_<data_suffix_U>.nc
296 FORMAT(a,'/',a,'_',i4.4,i2.2,'_',a,'.nc')                        ! <data_dir>/<data_prefix>_YYYYMM_<data_suffix_U>.nc
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
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_gridU,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_gridU,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
        CASE(196) 
          write(file_in_gridU,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U)
        CASE(197) 
          write(file_in_gridU,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_gridU,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_gridU,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
        CASE(292)
          write(file_in_gridU,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U) 
        CASE(293)
          write(file_in_gridU,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_gridU,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_gridU,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
        CASE(296) 
          write(file_in_gridU,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U)
        CASE(297) 
          write(file_in_gridU,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_gridU,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     inquire(file=file_in_gridU, exist=existfile)
     if ( existfile ) exit
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading U input dimensions in ', TRIM(file_in_gridU)
    status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)          ; call erreur(status,.TRUE.,"read first")

    status = NF90_INQ_DIMID(fidU,"time_counter",dimID_time) ; call erreur(status,.TRUE.,"inq_dimID_time")
    status = NF90_INQ_DIMID(fidU,"x",dimID_x)               ; call erreur(status,.TRUE.,"inq_dimID_x")
    status = NF90_INQ_DIMID(fidU,"y",dimID_y)               ; call erreur(status,.TRUE.,"inq_dimID_y")
    status = NF90_INQ_DIMID(fidU,"z",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depth",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"depthv",dimID_depthv)
    if (status .ne. 0) status = NF90_INQ_DIMID(fidU,"deptht",dimID_depthv)
    call erreur(status,.TRUE.,"inq_dimID_depthv")

    status = NF90_INQUIRE_DIMENSION(fidU,dimID_time,len=mtime)     ; call erreur(status,.TRUE.,"inq_dim_time")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_x,len=mlon)         ; call erreur(status,.TRUE.,"inq_dim_x")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_y,len=mlat)         ; call erreur(status,.TRUE.,"inq_dim_y")
    status = NF90_INQUIRE_DIMENSION(fidU,dimID_depthv,len=mdepthv) ; call erreur(status,.TRUE.,"inq_dim_depthv")

    ALLOCATE( nav_lon(mlon,mlat), nav_lat(mlon,mlat) )
    ALLOCATE( depthv(mdepthv) )

    status = NF90_INQ_VARID(fidU,"nav_lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lon",lon_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"longitude",lon_ID)
    call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_INQ_VARID(fidU,"nav_lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"lat",lat_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"latitude",lat_ID)
    call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidU,"depthv",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"depth",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"deptht",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"nav_lev",depth_ID)
    if ( status .ne. 0 ) status = NF90_INQ_VARID(fidU,"z",depth_ID)
    call erreur(status,.TRUE.,"inq_depthv_ID")
        
    status = NF90_GET_VAR(fidU,lon_ID,nav_lon)                   ; call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_VAR(fidU,lat_ID,nav_lat)                   ; call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidU,depth_ID,depthv)                  ; call erreur(status,.TRUE.,"getvar_depthv")

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
! 3- Process all gridU files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_gridU,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
          write(file_in_gridV,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_V)
        CASE(192)
          write(file_in_gridU,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U) 
          write(file_in_gridV,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_V) 
        CASE(193)
          write(file_in_gridU,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridV,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_gridU,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
          write(file_in_gridV,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_gridU,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
          write(file_in_gridV,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_V)
        CASE(196) 
          write(file_in_gridU,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U)
          write(file_in_gridV,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_V)
        CASE(197) 
          write(file_in_gridU,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridV,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_gridU,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
          write(file_in_gridV,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_gridU,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
          write(file_in_gridV,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_V)
        CASE(292)
          write(file_in_gridU,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U) 
          write(file_in_gridV,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_V) 
        CASE(293)
          write(file_in_gridU,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridV,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_gridU,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
          write(file_in_gridV,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_gridU,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_U)
          write(file_in_gridV,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_V)
        CASE(296) 
          write(file_in_gridU,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_U)
          write(file_in_gridV,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_V)
        CASE(297) 
          write(file_in_gridU,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
          write(file_in_gridV,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_gridU,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
          write(file_in_gridV,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_gridU, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          401 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV2d,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          501 FORMAT(a,'/BDY/bdyV_u3d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV3d,501) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          402 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV2d,402) TRIM(config_dir), kyear, kmonth, TRIM(config)
          502 FORMAT(a,'/BDY/bdyV_u3d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridV3d,502) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_gridV2d and file_bdy_gridV3d  >>>> stop'
          stop
        endif

        ALLOCATE( vomecrty_GLO(mlon,mlat,mdepthv,mtime)  )
        ALLOCATE( vozocrtx_GLO(mlon,mlat,mdepthv,mtime)  )
        ALLOCATE( vozocrtx_GLO_V(mlon,mlat,mdepthv,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input velocity :

        write(*,*) 'Reading U in ', TRIM(file_in_gridU)
        
        status = NF90_OPEN(TRIM(file_in_gridU),0,fidU)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidU,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidU,"vozocrtx",vozocrtx_ID)          ; call erreur(status,.TRUE.,"inq_vozocrtx_ID")
        
        status = NF90_GET_VAR(fidU,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidU,vozocrtx_ID,vozocrtx_GLO)          ; call erreur(status,.TRUE.,"getvar_vozocrtx")

        status = NF90_GET_ATT(fidU,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidU,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidU)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Read input velocity :

        write(*,*) 'Reading V in ', TRIM(file_in_gridV)
        
        status = NF90_OPEN(TRIM(file_in_gridV),0,fidV)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidV,"vomecrty",vomecrty_ID)          ; call erreur(status,.TRUE.,"inq_vomecrty_ID")
        
        status = NF90_GET_VAR(fidV,vomecrty_ID,vomecrty_GLO)          ; call erreur(status,.TRUE.,"getvar_vomecrty")

        status = NF90_CLOSE(fidV)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :

        write(*,*) 'Remove NaNs'
 
        do i=1,mlon
        do j=1,mlat
        do k=1,mdepthv
        do l=1,mtime
          if ( .not. abs(vozocrtx_GLO(i,j,k,l)) .lt. 100.0 ) then
            vozocrtx_GLO(i,j,k,l) = 0.0
          endif
          if ( .not. abs(vomecrty_GLO(i,j,k,l)) .lt. 100.0 ) then
            vomecrty_GLO(i,j,k,l) = 0.0
          endif
        enddo
        enddo
        enddo
        enddo

        !---------------------------------------
        ! Express velocities onto the other grid :

        write(*,*) 'Express V on grudU'

        do i=1,mlon
        do j=1,mlat

          ip1=MIN(i+1,mlon)
          im1=MAX(i-1,1)
          jp1=MIN(j+1,mlat)
          jm1=MAX(j-1,1)

          vozocrtx_GLO_V(i,j,:,:) = (   vozocrtx_GLO(i  ,j  ,:,:) * e2u_GLO(i  ,j  ) &
          &                           + vozocrtx_GLO(i  ,jp1,:,:) * e2u_GLO(i  ,jp1) &
          &                           + vozocrtx_GLO(im1,j  ,:,:) * e2u_GLO(im1,j  ) &
          &                           + vozocrtx_GLO(im1,jp1,:,:) * e2u_GLO(im1,jp1) ) / (4*e2v_GLO(i,j)) 

        enddo
        enddo

        !---------------------------------------
        ! Fill values on bdyU :
      
        write(*,*) 'Interpolating then rotating...'

        ALLOCATE( vomecrty_bdy(mxbv,1,mdepthv,mtime)  )

        do kbdy=1,mxbv

          do kz=1,mdepthv

            aNW = wgNWv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz ) &
                                                   * e3v_GLO  ( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz )
            aSW = wgSWv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz ) &
                                                   * e3v_GLO  ( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz )
            aNE = wgNEv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz ) &
                                                   * e3v_GLO  ( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz )
            aSE = wgSEv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz ) &
                                                   * e3v_GLO  ( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz )
             
            aa = aNW + aSW + aNE + aSE

            if ( aa .gt. eps .and. zjSv(nbiv(kbdy,1),nbjv(kbdy,1)) .gt. 1 ) then

              do kt=1,mtime

                vv = (   vomecrty_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aNW   &
                &      + vomecrty_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aSW   &
                &      + vomecrty_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aNE   &
                &      + vomecrty_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aSE ) &
                &    * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz) / aa

                uu = (   vozocrtx_GLO_V( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aNW   &
                &      + vozocrtx_GLO_V( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aSW   &
                &      + vozocrtx_GLO_V( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aNE   &
                &      + vozocrtx_GLO_V( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kz, kt ) * aSE ) &
                &    * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz) / aa

                vomecrty_bdy(kbdy,1,kz,kt) = - uu * sin(angXv(nbiv(kbdy,1),nbjv(kbdy,1))) + vv * cos(angYv(nbiv(kbdy,1),nbjv(kbdy,1)))

              enddo !- kt

            else
 
              vomecrty_bdy(kbdy,1,kz,:) = 0.e0
  
            endif

          enddo !- kz

          !- Removing barotropic component :
          do kt=1,mtime
            tmpdep = 0.e0
            tmpbtp = 0.e0
            do kz=1,mdepthv                   
               tmpbtp = tmpbtp + vomecrty_bdy(kbdy,1,kz,kt) * e3v_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz) * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz)
               tmpdep = tmpdep +                              e3v_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz) * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz)
            enddo
            if ( tmpdep .gt. eps ) tmpbtp = tmpbtp / tmpdep
            do kz=1,mdepthv
              vomecrty_bdy(kbdy,1,kz,kt) = ( vomecrty_bdy(kbdy,1,kz,kt) - tmpbtp ) * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),kz)
            enddo
          enddo

        enddo !- kbdy

        !--------------------------------------
        ! Write BDY netcdf file for baroclinic velocity :

        write(*,*) 'Creating ', TRIM(file_bdy_gridV3d)
        status = NF90_CREATE(TRIM(file_bdy_gridV3d),NF90_NOCLOBBER,fidN) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidN,"depthv",mdepthv,dimID_depthv)                    ; call erreur(status,.TRUE.,"def_dimID_depthv")
        status = NF90_DEF_DIM(fidN,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidN,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidN,"xbV",mxbV,dimID_xbV)                             ; call erreur(status,.TRUE.,"def_dimID_xbV")

        status = NF90_DEF_VAR(fidN,"vomecrty",NF90_FLOAT,(/dimID_xbV,dimID_yb,dimID_depthv,dimID_time_counter/),vomecrty_ID)
        call erreur(status,.TRUE.,"def_var_vomecrty_ID")
        status = NF90_DEF_VAR(fidN,"depthv",NF90_FLOAT,(/dimID_depthv/),depthv_ID)
        call erreur(status,.TRUE.,"def_var_depthv_ID")
        status = NF90_DEF_VAR(fidN,"nav_lat",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidN,"nav_lon",NF90_FLOAT,(/dimID_xbV,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidN,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidN,"nbrdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbrv_ID)
        call erreur(status,.TRUE.,"def_var_nbrv_ID")
        status = NF90_DEF_VAR(fidN,"nbjdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbjv_ID)
        call erreur(status,.TRUE.,"def_var_nbjv_ID")
        status = NF90_DEF_VAR(fidN,"nbidta",NF90_INT,(/dimID_xbV,dimID_yb/),nbiv_ID)
        call erreur(status,.TRUE.,"def_var_nbiv_ID")
       
        status = NF90_PUT_ATT(fidN,vomecrty_ID,"long_name","Baroclinic Velocity")   ; call erreur(status,.TRUE.,"put_att_vomecrty_ID")
        status = NF90_PUT_ATT(fidN,vomecrty_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vomecrty_ID")
        status = NF90_PUT_ATT(fidN,depthv_ID,"long_name","Vertical T levels")       ; call erreur(status,.TRUE.,"put_att_depthv_ID")
        status = NF90_PUT_ATT(fidN,depthv_ID,"units","m")                           ; call erreur(status,.TRUE.,"put_att_depthv_ID")
        status = NF90_PUT_ATT(fidN,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidN,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidN,nbrv_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidN,nbrv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidN,nbjv_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidN,nbjv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidN,nbiv_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        status = NF90_PUT_ATT(fidN,nbiv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"history","Created using extract_bdy_gridU.f90")
        status = NF90_PUT_ATT(fidN,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidN) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidN,vomecrty_ID,vomecrty_bdy) ; call erreur(status,.TRUE.,"var_vomecrty_ID")
        status = NF90_PUT_VAR(fidN,depthv_ID,depthv)         ; call erreur(status,.TRUE.,"var_depthv_ID")
        status = NF90_PUT_VAR(fidN,nav_lat_ID,gphiv_bdy)     ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidN,nav_lon_ID,glamv_bdy)     ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidN,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidN,nbrv_ID,nbrv)             ; call erreur(status,.TRUE.,"var_nbrv_ID")
        status = NF90_PUT_VAR(fidN,nbjv_ID,nbjv)             ; call erreur(status,.TRUE.,"var_nbjv_ID")
        status = NF90_PUT_VAR(fidN,nbiv_ID,nbiv)             ; call erreur(status,.TRUE.,"var_nbiv_ID")
        
        status = NF90_CLOSE(fidN) ; call erreur(status,.TRUE.,"close BDY file")
        
        !--       
        DEALLOCATE( vomecrty_GLO, vozocrtx_GLO, vozocrtx_GLO_V, time )
        DEALLOCATE( vomecrty_bdy )

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
