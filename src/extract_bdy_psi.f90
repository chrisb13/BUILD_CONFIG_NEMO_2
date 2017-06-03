program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to build netcdf coordinate file for BDY
!
! 0- Initialiartions
! 1- Read information on grids
! 2- Read input file dimensions in first existing file for specified time window
! 3- Process all PSI files over specified period
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
&                   data_suffix_bsf, data_suffix_ssh, data_suffix_ice,          &
&                   file_data_mask, file_data_zgr, file_data_hgr
CHARACTER(LEN=50)                    :: config
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_bsf, data_suffix_ice,  &
&                                       data_suffix_ssh, file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, myb, mxbt, glamt_ID, gphit_ID,   &
&                                       nbit_ID, nbjt_ID, nbrt_ID, mtime, dimID_x, dimID_y, vmask_REG_ID,       &
&                                       mlon, mlat, mdeptht, kday, kmonth, kyear, kbdy, nfmt, e1v_ID, e2u_ID,   &
&                                       kt, kz,     lon_ID, lat_ID, deptht_ID, fidM, nbiu_ID, nbju_ID, nbru_ID, &
&                                       sobarstf_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidT, dimID_xbu,  &
&                                       dimID_time_counter, dimID_deptht, time_ID, dimID_time, fidS, dimID_xbv, &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO, mxbu, mxbv,     &
&                                       depth_ID, ai, aj, bi, bj, kfmt, mx_REG, my_REG, umask_REG_ID, fidU,     &
&                                       tmask_GLO_ID, mx_tmp, my_tmp, fidMSKIN, fidZGR, e3t_GLO_ID, vobtcrtx_ID,&
&                                       mz_REG, fidMSKREG, dimID_z, mx_GLO, my_GLO, mz_GLO, vobtcrty_ID, fidV,  &
&                                       nbiv_ID, nbjv_ID, nbrv_ID, iREG, jREG, glamu_ID, glamv_ID, e3v_0_ID,    &
&                                       gphiu_ID, gphiv_ID, e3u_0_ID
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_PSI, file_bdy_PSI, file_bdy_gridU2d, &
&                                       file_in_coord_REG, command_str, file_in_mask_REG, file_bdy_gridV2d
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: tmask_GLO, umask_REG, vmask_REG
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt, nbiu, nbju, nbru, nbiv, nbjv, nbrv
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamu_bdy, gphiu_bdy, glamv_bdy, gphiv_bdy, nav_lon, nav_lat, e1v, e2u
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)  :: sobarstf_GLO, sobarstf_REG, vobtcrtx_bdy, vobtcrty_bdy, e3v_0, e3u_0
REAL*4,ALLOCATABLE,DIMENSION(:)      :: deptht
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
LOGICAL                              :: existfile

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEf_ID, wgNWf_ID, wgSWf_ID, wgSEf_ID, angYt_ID, angXt_ID,  &
&                                         zkUt_ID, ziWf_ID, ziEf_ID, zjNf_ID, zjSf_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWf, ziEf, zjNf, zjSf
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYt, angXt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEf, wgNWf, wgSWf, wgSEf

!-- Regional initial state
REAL*8                                 :: eps, aa

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

ALLOCATE(  tmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask_GLO); call erreur(status,.TRUE.,"getvar_tmask_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1c- Read BDY coordinates
!=================================================================================

write(*,*) 'Reading BDY coordinates in ', TRIM(file_coord)
status = NF90_OPEN(TRIM(file_coord),0,fidCOORD) ; call erreur(status,.TRUE.,"read bdy coordinate file") 

status = NF90_INQ_DIMID(fidCOORD,"yb",dimID_yb)   ; call erreur(status,.TRUE.,"inq_dimID_yb")
status = NF90_INQ_DIMID(fidCOORD,"xbt",dimID_xbt) ; call erreur(status,.TRUE.,"inq_dimID_xbt")
status = NF90_INQ_DIMID(fidCOORD,"xbu",dimID_xbu) ; call erreur(status,.TRUE.,"inq_dimID_xbu")
status = NF90_INQ_DIMID(fidCOORD,"xbv",dimID_xbv) ; call erreur(status,.TRUE.,"inq_dimID_xbv")

status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_yb,len=myb)   ; call erreur(status,.TRUE.,"inq_dim_yb")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbt,len=mxbt) ; call erreur(status,.TRUE.,"inq_dim_xbt")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbu,len=mxbu) ; call erreur(status,.TRUE.,"inq_dim_xbu")
status = NF90_INQUIRE_DIMENSION(fidCOORD,dimID_xbv,len=mxbv) ; call erreur(status,.TRUE.,"inq_dim_xbv")

ALLOCATE(  glamu_bdy(mxbu,myb)  ) 
ALLOCATE(  glamv_bdy(mxbv,myb)  ) 
ALLOCATE(  gphiu_bdy(mxbu,myb)  ) 
ALLOCATE(  gphiv_bdy(mxbv,myb)  ) 
ALLOCATE(  nbit(mxbt,myb)  ) 
ALLOCATE(  nbjt(mxbt,myb)  ) 
ALLOCATE(  nbrt(mxbt,myb)  ) 
ALLOCATE(  nbiu(mxbu,myb)  )
ALLOCATE(  nbju(mxbu,myb)  )
ALLOCATE(  nbru(mxbu,myb)  )
ALLOCATE(  nbiv(mxbv,myb)  )
ALLOCATE(  nbjv(mxbv,myb)  )
ALLOCATE(  nbrv(mxbv,myb)  )

status = NF90_INQ_VARID(fidCOORD,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidCOORD,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidCOORD,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiu",nbiu_ID)   ; call erreur(status,.TRUE.,"inq_nbiu_ID")
status = NF90_INQ_VARID(fidCOORD,"nbju",nbju_ID)   ; call erreur(status,.TRUE.,"inq_nbju_ID")
status = NF90_INQ_VARID(fidCOORD,"nbru",nbru_ID)   ; call erreur(status,.TRUE.,"inq_nbru_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiv",nbiv_ID)   ; call erreur(status,.TRUE.,"inq_nbiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjv",nbjv_ID)   ; call erreur(status,.TRUE.,"inq_nbjv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrv",nbrv_ID)   ; call erreur(status,.TRUE.,"inq_nbrv_ID")

status = NF90_GET_VAR(fidCOORD,glamu_ID,glamu_bdy) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidCOORD,glamv_ID,glamv_bdy) ; call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidCOORD,gphiu_ID,gphiu_bdy) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidCOORD,gphiv_ID,gphiv_bdy) ; call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidCOORD,nbit_ID,nbit)       ; call erreur(status,.TRUE.,"getvar_nbit")
status = NF90_GET_VAR(fidCOORD,nbjt_ID,nbjt)       ; call erreur(status,.TRUE.,"getvar_nbjt")
status = NF90_GET_VAR(fidCOORD,nbrt_ID,nbrt)       ; call erreur(status,.TRUE.,"getvar_nbrt")
status = NF90_GET_VAR(fidCOORD,nbiu_ID,nbiu)       ; call erreur(status,.TRUE.,"getvar_nbiu")
status = NF90_GET_VAR(fidCOORD,nbju_ID,nbju)       ; call erreur(status,.TRUE.,"getvar_nbju")
status = NF90_GET_VAR(fidCOORD,nbru_ID,nbru)       ; call erreur(status,.TRUE.,"getvar_nbru")
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

ALLOCATE(  umask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  vmask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  e3v_0    (mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  e3u_0    (mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  e2u      (mx_REG,my_REG       )  ) 
ALLOCATE(  e1v      (mx_REG,my_REG       )  ) 

status = NF90_INQ_VARID(fidMSKREG,"umask",umask_REG_ID); call erreur(status,.TRUE.,"inq_umask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"vmask",vmask_REG_ID); call erreur(status,.TRUE.,"inq_vmask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"e3v_0",e3v_0_ID);     call erreur(status,.TRUE.,"inq_e3v_0_ID")
status = NF90_INQ_VARID(fidMSKREG,"e3u_0",e3u_0_ID);     call erreur(status,.TRUE.,"inq_e3u_0_ID")
status = NF90_INQ_VARID(fidMSKREG,"e2u",e2u_ID);         call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidMSKREG,"e1v",e1v_ID);         call erreur(status,.TRUE.,"inq_e1v_ID")

status = NF90_GET_VAR(fidMSKREG,umask_REG_ID,umask_REG); call erreur(status,.TRUE.,"getvar_umask_REG")
status = NF90_GET_VAR(fidMSKREG,vmask_REG_ID,vmask_REG); call erreur(status,.TRUE.,"getvar_vmask_REG")
status = NF90_GET_VAR(fidMSKREG,e3v_0_ID,e3v_0);         call erreur(status,.TRUE.,"getvar_e3v_0")
status = NF90_GET_VAR(fidMSKREG,e3u_0_ID,e3u_0);         call erreur(status,.TRUE.,"getvar_e3u_0")
status = NF90_GET_VAR(fidMSKREG,e2u_ID,e2u);             call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidMSKREG,e1v_ID,e1v);             call erreur(status,.TRUE.,"getvar_e1v")

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

ALLOCATE(  wgNEf(mx_REG,my_REG), wgNWf(mx_REG,my_REG), wgSWf(mx_REG,my_REG), wgSEf(mx_REG,my_REG) )
ALLOCATE(  ziWf(mx_REG,my_REG), ziEf(mx_REG,my_REG), zjNf(mx_REG,my_REG), zjSf(mx_REG,my_REG)  ) 
            
status = NF90_INQ_VARID(fidcoeff,"wgNEf",wgNEf_ID) ; call erreur(status,.TRUE.,"inq_wgNEf_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWf",wgNWf_ID) ; call erreur(status,.TRUE.,"inq_wgNWf_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWf",wgSWf_ID) ; call erreur(status,.TRUE.,"inq_wgSWf_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEf",wgSEf_ID) ; call erreur(status,.TRUE.,"inq_wgSEf_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNf",zjNf_ID)   ; call erreur(status,.TRUE.,"inq_zjNf_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSf",zjSf_ID)   ; call erreur(status,.TRUE.,"inq_zjSf_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEf",ziEf_ID)   ; call erreur(status,.TRUE.,"inq_ziEf_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWf",ziWf_ID)   ; call erreur(status,.TRUE.,"inq_ziWf_ID")
               
status = NF90_GET_VAR(fidcoeff,wgNEf_ID,wgNEf)     ; call erreur(status,.TRUE.,"getvar_wgNEf")
status = NF90_GET_VAR(fidcoeff,wgNWf_ID,wgNWf)     ; call erreur(status,.TRUE.,"getvar_wgNWf")
status = NF90_GET_VAR(fidcoeff,wgSWf_ID,wgSWf)     ; call erreur(status,.TRUE.,"getvar_wgSWf")
status = NF90_GET_VAR(fidcoeff,wgSEf_ID,wgSEf)     ; call erreur(status,.TRUE.,"getvar_wgSEf")
status = NF90_GET_VAR(fidcoeff,zjNf_ID,zjNf)       ; call erreur(status,.TRUE.,"getvar_zjNf")
status = NF90_GET_VAR(fidcoeff,zjSf_ID,zjSf)       ; call erreur(status,.TRUE.,"getvar_zjSf")
status = NF90_GET_VAR(fidcoeff,ziEf_ID,ziEf)       ; call erreur(status,.TRUE.,"getvar_ziEf")
status = NF90_GET_VAR(fidcoeff,ziWf_ID,ziWf)       ; call erreur(status,.TRUE.,"getvar_ziWf")
                            
status = NF90_CLOSE(fidcoeff) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 4- Read input file dimensions in first existing file for specified time window
!=================================================================================

!- accepted input format :
191 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')  ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD_<data_suffix_bsf>.nc
192 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')           ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_<data_suffix_bsf>.nc
193 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')        ! <data_dir>/YYYY/<data_prefix>_YYYY_MM_DD.nc
194 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,'_',i2.2,'.nc')                 ! <data_dir>/YYYY/<data_prefix>_YYYY_MM.nc
195 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')           ! <data_dir>/<data_prefix>_YYYY_MM_DD_<data_suffix_bsf>.nc
196 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',a,'.nc')                    ! <data_dir>/<data_prefix>_YYYY_MM_<data_suffix_bsf>.nc
197 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'_',i2.2,'.nc')                 ! <data_dir>/<data_prefix>_YYYY_MM_DD.nc
198 FORMAT(a,'/',a,'_',i4.4,'_',i2.2,'.nc')                          ! <data_dir>/<data_prefix>_YYYY_MM.nc
291 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')          ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD_<data_suffix_bsf>.nc
292 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'_',a,'.nc')               ! <data_dir>/YYYY/<data_prefix>_YYYYMM_<data_suffix_bsf>.nc
293 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,i2.2,'.nc')                ! <data_dir>/YYYY/<data_prefix>_YYYYMMDD.nc
294 FORMAT(a,'/',i4.4,'/',a,'_',i4.4,i2.2,'.nc')                     ! <data_dir>/YYYY/<data_prefix>_YYYYMM.nc
295 FORMAT(a,'/',a,'_',i4.4,i2.2,i2.2,'_',a,'.nc')                   ! <data_dir>/<data_prefix>_YYYYMMDD_<data_suffix_bsf>.nc
296 FORMAT(a,'/',a,'_',i4.4,i2.2,'_',a,'.nc')                        ! <data_dir>/<data_prefix>_YYYYMM_<data_suffix_bsf>.nc
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
          write(file_in_PSI,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(192)
          write(file_in_PSI,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf) 
        CASE(193)
          write(file_in_PSI,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_PSI,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_PSI,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(196) 
          write(file_in_PSI,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf)
        CASE(197) 
          write(file_in_PSI,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_PSI,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_PSI,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(292)
          write(file_in_PSI,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf) 
        CASE(293)
          write(file_in_PSI,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_PSI,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_PSI,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(296) 
          write(file_in_PSI,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf)
        CASE(297) 
          write(file_in_PSI,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_PSI,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     inquire(file=file_in_PSI, exist=existfile)
     if ( existfile ) exit
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_PSI)
    status = NF90_OPEN(TRIM(file_in_PSI),0,fidT)          ; call erreur(status,.TRUE.,"read first")

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

    write(*,*) 'No BSF file found for first month of ', kyear
    write(*,*) '      >>>>>>>>>>>>>>>>>> stop !!'
    stop

  ENDIF

ENDDO

!--

write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a,'/BDY')
CALL system(TRIM(command_str))

!=================================================================================
! 3- Process all PSI files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_PSI,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(192)
          write(file_in_PSI,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf) 
        CASE(193)
          write(file_in_PSI,193) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(194) 
          write(file_in_PSI,194) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(195) 
          write(file_in_PSI,195) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(196) 
          write(file_in_PSI,196) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf)
        CASE(197) 
          write(file_in_PSI,197) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(198)
          write(file_in_PSI,198) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE(291)
          write(file_in_PSI,291) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(292)
          write(file_in_PSI,292) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf) 
        CASE(293)
          write(file_in_PSI,293) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth, kday
        CASE(294) 
          write(file_in_PSI,294) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, kmonth
        CASE(295) 
          write(file_in_PSI,295) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday, TRIM(data_suffix_bsf)
        CASE(296) 
          write(file_in_PSI,296) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, TRIM(data_suffix_bsf)
        CASE(297) 
          write(file_in_PSI,297) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth, kday
        CASE(298)
          write(file_in_PSI,298) TRIM(data_dir), TRIM(data_prefix), kyear, kmonth
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_PSI, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 .or. nfmt .eq. 195 .or. nfmt .eq. 197 &
        &   .or. nfmt .eq. 291 .or. nfmt .eq. 293 .or. nfmt .eq. 295 .or. nfmt .eq. 297 ) then
          401 FORMAT(a,'/BDY/bdyU_u2d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          402 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU2d,401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          write(file_bdy_gridV2d,402) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 .or. nfmt .eq. 196 .or. nfmt .eq. 198 &
        &   .or. nfmt .eq. 292 .or. nfmt .eq. 294 .or. nfmt .eq. 296 .or. nfmt .eq. 298 ) then
          403 FORMAT(a,'/BDY/bdyU_u2d_',i4.4,'_',i2.2,'_',a,'.nc')
          404 FORMAT(a,'/BDY/bdyV_u2d_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_gridU2d,403) TRIM(config_dir), kyear, kmonth, TRIM(config)
          write(file_bdy_gridV2d,404) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_PSI  >>>> stop'
          stop
        endif

        ALLOCATE( sobarstf_GLO(mlon,mlat,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! Read input psi :

        write(*,*) 'Reading psi in ', TRIM(file_in_PSI)
        
        status = NF90_OPEN(TRIM(file_in_PSI),0,fidT)                ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidT,"time_counter",time_ID)          ; call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidT,"sobarstf",sobarstf_ID)          ; call erreur(status,.TRUE.,"inq_sobarstf_ID")
        
        status = NF90_GET_VAR(fidT,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidT,sobarstf_ID,sobarstf_GLO)          ; call erreur(status,.TRUE.,"getvar_sobarstf")

        status = NF90_GET_ATT(fidT,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidT,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! Remove possible NaNs :
 
        do i=1,mlon
        do j=1,mlat
        do l=1,mtime
          if ( .not. abs(sobarstf_GLO(i,j,l)) .lt. 1.e10 ) then
            !sobarstf_GLO(i,j,l) = 0.0
            write(*,*) '~!@#$%^* ERROR: the BSF is assumed to have values everywhere >>>>>> stop !'
            stop
          endif
        enddo
        enddo
        enddo

        !------------------------------------------------
        ! Interpolate BSF values on the regional domain :
      
        ALLOCATE( sobarstf_REG(mx_REG,my_REG,mtime)  )
 
        do iREG=1,mx_REG
        do jREG=1,my_REG 

          aa = wgNWf(iREG,jREG) + wgSWf(iREG,jREG) + wgNEf(iREG,jREG) + wgSEf(iREG,jREG)

          if ( aa .gt. eps .and. zjSf(iREG,jREG) .gt. 1 ) then

            do kt=1,mtime
               sobarstf_REG (iREG,jREG,kt) =  (   sobarstf_GLO( ziWf(iREG,jREG), zjNf(iREG,jREG), kt ) * wgNWf(iREG,jREG)   &
               &                                + sobarstf_GLO( ziWf(iREG,jREG), zjSf(iREG,jREG), kt ) * wgSWf(iREG,jREG)   &
               &                                + sobarstf_GLO( ziEf(iREG,jREG), zjNf(iREG,jREG), kt ) * wgNEf(iREG,jREG)   &
               &                                + sobarstf_GLO( ziEf(iREG,jREG), zjSf(iREG,jREG), kt ) * wgSEf(iREG,jREG) ) / aa
            enddo

          else
 
            sobarstf_REG (iREG,jREG,:) = 0.0  !! should not be used anyway... 
  
          endif

        enddo
        enddo

        !--------------------------------------
        ! Extract bdyU_u2d from the BSF :

        ALLOCATE( vobtcrtx_bdy(mxbu,1,mtime)  )

        vobtcrtx_bdy(:,:,:) = 0.e0

        do kbdy=1,mxbu

         aa = e2u(nbiu(kbdy,1),nbju(kbdy,1)) * sum( e3u_0(nbiu(kbdy,1),nbju(kbdy,1),:)*umask_REG(nbiu(kbdy,1),nbju(kbdy,1),:) )

         ! U[i,j] = - ( PSIF[i,j] - PSIF[i,j-1] ) / aa
         if ( nbju(kbdy,1) .gt. 1 .and. aa .gt. eps ) then
           vobtcrtx_bdy(kbdy,1,:) = ( sobarstf_REG (nbiu(kbdy,1),nbju(kbdy,1)-1,:) - sobarstf_REG (nbiu(kbdy,1),nbju(kbdy,1),:) ) / aa
         endif
 
        enddo

        !--------------------------------------
        ! Extract bdyV_u2d from the BSF :

        ALLOCATE( vobtcrty_bdy(mxbv,1,mtime)  )

        vobtcrty_bdy(:,:,:) = 0.e0

        do kbdy=1,mxbv

         aa = e1v(nbiv(kbdy,1),nbjv(kbdy,1)) * sum( e3v_0(nbiv(kbdy,1),nbjv(kbdy,1),:)*vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),:) )

         ! V[i,j] = ( PSIF[i,j] - PSIF[i-1,j] ) / aa
         if ( nbiv(kbdy,1) .gt. 1 .and. aa .gt. eps ) then
           vobtcrty_bdy(kbdy,1,:) = ( sobarstf_REG (nbiv(kbdy,1),nbjv(kbdy,1),:) - sobarstf_REG (nbiv(kbdy,1)-1,nbjv(kbdy,1),:) ) / aa
         endif

        enddo

        !------------------------------------------------
        ! Write BDY netcdf file for U barotropic velocity

        write(*,*) 'Creating ', TRIM(file_bdy_gridU2d)
        status = NF90_CREATE(TRIM(file_bdy_gridU2d),NF90_NOCLOBBER,fidU) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidU,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidU,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidU,"xbu",mxbu,dimID_xbu)                             ; call erreur(status,.TRUE.,"def_dimID_xbu")

        status = NF90_DEF_VAR(fidU,"vobtcrtx",NF90_FLOAT,(/dimID_xbu,dimID_yb,dimID_time_counter/),vobtcrtx_ID)
        call erreur(status,.TRUE.,"def_var_vobtcrtx_ID")
        status = NF90_DEF_VAR(fidU,"nav_lat",NF90_FLOAT,(/dimID_xbu,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidU,"nav_lon",NF90_FLOAT,(/dimID_xbu,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidU,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidU,"nbrdta",NF90_INT,(/dimID_xbu,dimID_yb/),nbru_ID)
        call erreur(status,.TRUE.,"def_var_nbru_ID")
        status = NF90_DEF_VAR(fidU,"nbjdta",NF90_INT,(/dimID_xbu,dimID_yb/),nbju_ID)
        call erreur(status,.TRUE.,"def_var_nbju_ID")
        status = NF90_DEF_VAR(fidU,"nbidta",NF90_INT,(/dimID_xbu,dimID_yb/),nbiu_ID)
        call erreur(status,.TRUE.,"def_var_nbiu_ID")
       
        status = NF90_PUT_ATT(fidU,vobtcrtx_ID,"long_name","Barotropic velocity")   ; call erreur(status,.TRUE.,"put_att_vobtcrtx_ID")
        status = NF90_PUT_ATT(fidU,vobtcrtx_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vobtcrtx_ID")
        status = NF90_PUT_ATT(fidU,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidU,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidU,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidU,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidU,nbru_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidU,nbru_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidU,nbju_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidU,nbju_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidU,nbiu_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        status = NF90_PUT_ATT(fidU,nbiu_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        
        status = NF90_PUT_ATT(fidU,NF90_GLOBAL,"history","Created using extract_bdy_psi.f90")
        status = NF90_PUT_ATT(fidU,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidU) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidU,vobtcrtx_ID,vobtcrtx_bdy) ; call erreur(status,.TRUE.,"var_vobtcrtx_ID")
        status = NF90_PUT_VAR(fidU,nav_lat_ID,gphiu_bdy)     ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidU,nav_lon_ID,glamu_bdy)     ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidU,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidU,nbru_ID,nbru)             ; call erreur(status,.TRUE.,"var_nbru_ID")
        status = NF90_PUT_VAR(fidU,nbju_ID,nbju)             ; call erreur(status,.TRUE.,"var_nbju_ID")
        status = NF90_PUT_VAR(fidU,nbiu_ID,nbiu)             ; call erreur(status,.TRUE.,"var_nbiu_ID")
        
        status = NF90_CLOSE(fidU) ; call erreur(status,.TRUE.,"close BDY file")

        !------------------------------------------------
        ! Write BDY netcdf file for V barotropic velocity

        write(*,*) 'Creating ', TRIM(file_bdy_gridV2d)
        status = NF90_CREATE(TRIM(file_bdy_gridV2d),NF90_NOCLOBBER,fidV) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidV,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidV,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidV,"xbv",mxbv,dimID_xbv)                             ; call erreur(status,.TRUE.,"def_dimID_xbv")

        status = NF90_DEF_VAR(fidV,"vobtcrty",NF90_FLOAT,(/dimID_xbv,dimID_yb,dimID_time_counter/),vobtcrty_ID)
        call erreur(status,.TRUE.,"def_var_vobtcrty_ID")
        status = NF90_DEF_VAR(fidV,"nav_lat",NF90_FLOAT,(/dimID_xbv,dimID_yb/),nav_lat_ID)
        call erreur(status,.TRUE.,"def_var_nav_lat_ID")
        status = NF90_DEF_VAR(fidV,"nav_lon",NF90_FLOAT,(/dimID_xbv,dimID_yb/),nav_lon_ID)
        call erreur(status,.TRUE.,"def_var_nav_lon_ID")
        status = NF90_DEF_VAR(fidV,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidV,"nbrdta",NF90_INT,(/dimID_xbv,dimID_yb/),nbrv_ID)
        call erreur(status,.TRUE.,"def_var_nbrv_ID")
        status = NF90_DEF_VAR(fidV,"nbjdta",NF90_INT,(/dimID_xbv,dimID_yb/),nbjv_ID)
        call erreur(status,.TRUE.,"def_var_nbjv_ID")
        status = NF90_DEF_VAR(fidV,"nbidta",NF90_INT,(/dimID_xbv,dimID_yb/),nbiv_ID)
        call erreur(status,.TRUE.,"def_var_nbiv_ID")
       
        status = NF90_PUT_ATT(fidV,vobtcrty_ID,"long_name","Barotropic velocity")   ; call erreur(status,.TRUE.,"put_att_vobtcrty_ID")
        status = NF90_PUT_ATT(fidV,vobtcrty_ID,"units","m/s")                       ; call erreur(status,.TRUE.,"put_att_vobtcrty_ID")
        status = NF90_PUT_ATT(fidV,nav_lat_ID,"units","degrees_north")              ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidV,nav_lon_ID,"units","degrees_east")               ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidV,time_counter_ID,"units",TRIM(time_units))        ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidV,time_counter_ID,"calendar",TRIM(calendar))       ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidV,nbrv_ID,"long_name","bdy discrete distance")     ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidV,nbrv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidV,nbjv_ID,"long_name","bdy j index")               ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidV,nbjv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidV,nbiv_ID,"long_name","bdy i index")               ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        status = NF90_PUT_ATT(fidV,nbiv_ID,"units","unitless")                      ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        
        status = NF90_PUT_ATT(fidV,NF90_GLOBAL,"history","Created using extract_bdy_psi.f90")
        status = NF90_PUT_ATT(fidV,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidV) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidV,vobtcrty_ID,vobtcrty_bdy) ; call erreur(status,.TRUE.,"var_vobtcrty_ID")
        status = NF90_PUT_VAR(fidV,nav_lat_ID,gphiv_bdy)     ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidV,nav_lon_ID,glamv_bdy)     ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidV,time_counter_ID,time)     ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidV,nbrv_ID,nbrv)             ; call erreur(status,.TRUE.,"var_nbrv_ID")
        status = NF90_PUT_VAR(fidV,nbjv_ID,nbjv)             ; call erreur(status,.TRUE.,"var_nbjv_ID")
        status = NF90_PUT_VAR(fidV,nbiv_ID,nbiv)             ; call erreur(status,.TRUE.,"var_nbiv_ID")
        
        status = NF90_CLOSE(fidV) ; call erreur(status,.TRUE.,"close BDY file")

        
        !--       
        DEALLOCATE( sobarstf_GLO, time )
        DEALLOCATE( sobarstf_REG, vobtcrtx_bdy, vobtcrty_bdy )

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