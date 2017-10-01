program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, March 2015
!
! Used to extract SSH along the BDY
!
! 0- Initialiartions
! 1a- Read GLOBAL mask
! 1b- Read GLOBAL horizontal grid metrics
! 1c- Read BDY coordinates
! 2- Read REGIONAL mask
! 3- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
! 4- Read input file dimensions in first existing file for specified time window
! 5- Process all SSH files over specified period
! 5a- Read input (GLOBAL) ssh
! 5b- Calculate SSH on U and V grids
! 5c- Remove possible NaNs
! 5d- Calculate SSH values on bdyT
! 5e- Calculate SSH values on bdyU
! 5f- Calculate SSH values on bdyV
! 5g- Write new BDY netcdf file containing SSH at T points
! 5h- Write new BDY netcdf file containing SSH at U points
! 5i- Write new BDY netcdf file containing SSH at V points
!
! history : - Feb. 2017: version with namelist (nj)
!           - Sep. 2017: save SSH on grids U and V as well (useful to derive the
!                        BSF in extract_bdy_psi )
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
CHARACTER(LEN=150)                   :: config_dir, data_dir, data_prefix, data_suffix_T, data_suffix_S, &
&                                       data_suffix_U, data_suffix_V, data_suffix_ssh, data_suffix_ice,  &
&                                       data_suffix_bsf, file_data_mask, file_data_zgr, file_data_hgr
INTEGER                              :: nn_yeari, nn_yearf, nn_bdy_eosmatch

INTEGER                              :: fidCOORD, status, dimID_yb, dimID_xbt, myb, mxbt, glamt_ID, gphit_ID, &
&                                       e1t_ID, e2t_ID, nbit_ID, nbjt_ID, nbrt_ID, mtime, dimID_x, dimID_y,   &
&                                       mlon, mlat, mdeptht, kday, kmonth, kyear, kbdy, nfmt, e1u_ID, e2u_ID, &
&                                       kt, kz,     lon_ID, lat_ID, deptht_ID, fidM, e1v_ID, e2v_ID, mx_GLO,  &
&                                       ssh_ID, time_counter_ID, nav_lon_ID, nav_lat_ID, fidT, umask_GLO_ID,&
&                                       dimID_time_counter, dimID_deptht, time_ID, dimID_time, fidS, my_GLO,  &
&                                       i, j, k, l, fidC, imin_ORCA12, jmin_ORCA12, iGLO, jGLO, vmask_GLO_ID, &
&                                       depth_ID, ai, aj, bi, bj, kfmt, mx_REG, my_REG, tmask_REG_ID, mz_GLO, &
&                                       tmask_GLO_ID, mx_tmp, my_tmp, fidMSKIN, fidZGR, e3t_GLO_ID, mz_REG,   &
&                                       fidMSKREG, dimID_z, mxbu, mxbv, fidHGR, dimID_xbu, dimID_xbv, nbiu_ID,&
&                                       nbiv_ID, nbju_ID, nbjv_ID, nbru_ID, nbrv_ID, umask_REG_ID,            &
&                                       vmask_REG_ID
CHARACTER(LEN=100)                   :: calendar, time_units
CHARACTER(LEN=150)                   :: file_coord, file_in_SSH, file_bdy_SSH, file_bdy_SSH_U, file_bdy_SSH_V,&
&                                       file_in_coord_REG, command_str, file_in_mask_REG
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: tmask_GLO, tmask_REG, umask_GLO, umask_REG, vmask_GLO, vmask_REG
INTEGER*4,ALLOCATABLE,DIMENSION(:)   :: list_fmt
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: nbit, nbjt, nbrt, nbiu, nbju, nbru, nbiv, nbjv, nbrv
REAL*4,ALLOCATABLE,DIMENSION(:,:)    :: glamt_bdy, gphit_bdy, e1t, e1u, e1v, areaT, areaU, areaV, nav_lon, nav_lat
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)  :: ssh_T_GLO, ssh_T_bdy, ssh_U_GLO, ssh_U_bdy, ssh_V_GLO, ssh_V_bdy
REAL*4,ALLOCATABLE,DIMENSION(:)      :: deptht
REAL*8,ALLOCATABLE,DIMENSION(:)      :: time
LOGICAL                              :: existfile


!-- interpolation parameters
INTEGER                                :: fidcoeff, &
&                                         wgNEt_ID, wgNWt_ID, wgSWt_ID, wgSEt_ID,      &
&                                         zkUt_ID, ziWt_ID, ziEt_ID, zjNt_ID, zjSt_ID, &
&                                         wgNEu_ID, wgNWu_ID, wgSWu_ID, wgSEu_ID,      &
&                                         zkUu_ID, ziWu_ID, ziEu_ID, zjNu_ID, zjSu_ID, &
&                                         wgNEv_ID, wgNWv_ID, wgSWv_ID, wgSEv_ID,      &
&                                         zkUv_ID, ziWv_ID, ziEv_ID, zjNv_ID, zjSv_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWt, ziEt, zjNt, zjSt, &
&                                         ziWu, ziEu, zjNu, zjSu, &
&                                         ziWv, ziEv, zjNv, zjSv
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEt, wgNWt, wgSWt, wgSEt, &
&                                         wgNEu, wgNWu, wgSWu, wgSEu, &
&                                         wgNEv, wgNWv, wgSWv, wgSEv

!-- Regional initial state
REAL*8                                 :: eps, aa, aSW, aNW, aSE, aNE

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

ALLOCATE(  tmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 
ALLOCATE(  umask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 
ALLOCATE(  vmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")
status = NF90_INQ_VARID(fidMSKIN,"umask",umask_GLO_ID); call erreur(status,.TRUE.,"inq_umask_GLO_ID")
status = NF90_INQ_VARID(fidMSKIN,"vmask",vmask_GLO_ID); call erreur(status,.TRUE.,"inq_vmask_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask_GLO); call erreur(status,.TRUE.,"getvar_tmask_GLO")
status = NF90_GET_VAR(fidMSKIN,umask_GLO_ID,umask_GLO); call erreur(status,.TRUE.,"getvar_umask_GLO")
status = NF90_GET_VAR(fidMSKIN,vmask_GLO_ID,vmask_GLO); call erreur(status,.TRUE.,"getvar_vmask_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1b- Read GLOBAL horizontal grid metrics :
!=================================================================================

write(*,*) 'Reading ', TRIM(file_data_hgr)

status = NF90_OPEN(TRIM(file_data_hgr),0,fidHGR); call erreur(status,.TRUE.,"read hgr GLO") 

ALLOCATE(  e1t  (mx_GLO,my_GLO)  )
ALLOCATE(  areaT(mx_GLO,my_GLO)  )
ALLOCATE(  e1u  (mx_GLO,my_GLO)  )
ALLOCATE(  areaU(mx_GLO,my_GLO)  )
ALLOCATE(  e1v  (mx_GLO,my_GLO)  )
ALLOCATE(  areaV(mx_GLO,my_GLO)  )

status = NF90_INQ_VARID(fidHGR,"e1t",e1t_ID)     ; call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidHGR,"e2t",e2t_ID)     ; call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidHGR,"e1u",e1u_ID)     ; call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidHGR,"e2u",e2u_ID)     ; call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidHGR,"e1v",e1v_ID)     ; call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidHGR,"e2v",e2v_ID)     ; call erreur(status,.TRUE.,"inq_e2v_ID")

status = NF90_GET_VAR(fidHGR,e1t_ID,e1t)         ; call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidHGR,e2t_ID,areaT)       ; call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidHGR,e1u_ID,e1u)         ; call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidHGR,e2u_ID,areaU)       ; call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidHGR,e1v_ID,e1v)         ; call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidHGR,e2v_ID,areaV)       ; call erreur(status,.TRUE.,"getvar_e2v")

status = NF90_CLOSE(fidHGR); call erreur(status,.TRUE.,"end read hgr GLO")

areaT=e1t*areaT
areaU=e1u*areaU
areaV=e1v*areaV

DEALLOCATE( e1t, e1u, e1v )

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

ALLOCATE(  glamt_bdy(mxbt,myb)  ) 
ALLOCATE(  gphit_bdy(mxbt,myb)  ) 
ALLOCATE(  nbit (mxbt,myb)  ) 
ALLOCATE(  nbjt (mxbt,myb)  ) 
ALLOCATE(  nbrt (mxbt,myb)  ) 
ALLOCATE(  nbiu (mxbu,myb)  )
ALLOCATE(  nbju (mxbu,myb)  )
ALLOCATE(  nbru (mxbu,myb)  )
ALLOCATE(  nbiv (mxbv,myb)  )
ALLOCATE(  nbjv (mxbv,myb)  )
ALLOCATE(  nbrv (mxbv,myb)  )

status = NF90_INQ_VARID(fidCOORD,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidCOORD,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbit",nbit_ID)   ; call erreur(status,.TRUE.,"inq_nbit_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjt",nbjt_ID)   ; call erreur(status,.TRUE.,"inq_nbjt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrt",nbrt_ID)   ; call erreur(status,.TRUE.,"inq_nbrt_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiu",nbiu_ID)   ; call erreur(status,.TRUE.,"inq_nbiu_ID")
status = NF90_INQ_VARID(fidCOORD,"nbju",nbju_ID)   ; call erreur(status,.TRUE.,"inq_nbju_ID")
status = NF90_INQ_VARID(fidCOORD,"nbru",nbru_ID)   ; call erreur(status,.TRUE.,"inq_nbru_ID")
status = NF90_INQ_VARID(fidCOORD,"nbiv",nbiv_ID)   ; call erreur(status,.TRUE.,"inq_nbiv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbjv",nbjv_ID)   ; call erreur(status,.TRUE.,"inq_nbjv_ID")
status = NF90_INQ_VARID(fidCOORD,"nbrv",nbrv_ID)   ; call erreur(status,.TRUE.,"inq_nbrv_ID")

status = NF90_GET_VAR(fidCOORD,glamt_ID,glamt_bdy) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidCOORD,gphit_ID,gphit_bdy) ; call erreur(status,.TRUE.,"getvar_gphit")
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

ALLOCATE(  tmask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  umask_REG(mx_REG,my_REG,mz_REG)  ) 
ALLOCATE(  vmask_REG(mx_REG,my_REG,mz_REG)  ) 

status = NF90_INQ_VARID(fidMSKREG,"tmask",tmask_REG_ID); call erreur(status,.TRUE.,"inq_tmask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"umask",umask_REG_ID); call erreur(status,.TRUE.,"inq_umask_REG_ID")
status = NF90_INQ_VARID(fidMSKREG,"vmask",vmask_REG_ID); call erreur(status,.TRUE.,"inq_vmask_REG_ID")

status = NF90_GET_VAR(fidMSKREG,tmask_REG_ID,tmask_REG); call erreur(status,.TRUE.,"getvar_tmask_REG")
status = NF90_GET_VAR(fidMSKREG,umask_REG_ID,umask_REG); call erreur(status,.TRUE.,"getvar_umask_REG")
status = NF90_GET_VAR(fidMSKREG,vmask_REG_ID,vmask_REG); call erreur(status,.TRUE.,"getvar_vmask_REG")

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

ALLOCATE(  wgNEt(mx_REG,my_REG), wgNWt(mx_REG,my_REG), wgSWt(mx_REG,my_REG), wgSEt(mx_REG,my_REG)  )
ALLOCATE(  wgNEu(mx_REG,my_REG), wgNWu(mx_REG,my_REG), wgSWu(mx_REG,my_REG), wgSEu(mx_REG,my_REG)  )
ALLOCATE(  wgNEv(mx_REG,my_REG), wgNWv(mx_REG,my_REG), wgSWv(mx_REG,my_REG), wgSEv(mx_REG,my_REG)  )
ALLOCATE(   ziWt(mx_REG,my_REG),  ziEt(mx_REG,my_REG),  zjNt(mx_REG,my_REG),  zjSt(mx_REG,my_REG)  ) 
ALLOCATE(   ziWu(mx_REG,my_REG),  ziEu(mx_REG,my_REG),  zjNu(mx_REG,my_REG),  zjSu(mx_REG,my_REG)  ) 
ALLOCATE(   ziWv(mx_REG,my_REG),  ziEv(mx_REG,my_REG),  zjNv(mx_REG,my_REG),  zjSv(mx_REG,my_REG)  ) 
            
status = NF90_INQ_VARID(fidcoeff,"wgNEt",wgNEt_ID) ; call erreur(status,.TRUE.,"inq_wgNEt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWt",wgNWt_ID) ; call erreur(status,.TRUE.,"inq_wgNWt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWt",wgSWt_ID) ; call erreur(status,.TRUE.,"inq_wgSWt_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEt",wgSEt_ID) ; call erreur(status,.TRUE.,"inq_wgSEt_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNt",zjNt_ID)   ; call erreur(status,.TRUE.,"inq_zjNt_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSt",zjSt_ID)   ; call erreur(status,.TRUE.,"inq_zjSt_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEt",ziEt_ID)   ; call erreur(status,.TRUE.,"inq_ziEt_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWt",ziWt_ID)   ; call erreur(status,.TRUE.,"inq_ziWt_ID")
!-     
status = NF90_INQ_VARID(fidcoeff,"wgNEu",wgNEu_ID) ; call erreur(status,.TRUE.,"inq_wgNEu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWu",wgNWu_ID) ; call erreur(status,.TRUE.,"inq_wgNWu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWu",wgSWu_ID) ; call erreur(status,.TRUE.,"inq_wgSWu_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEu",wgSEu_ID) ; call erreur(status,.TRUE.,"inq_wgSEu_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNu",zjNu_ID)   ; call erreur(status,.TRUE.,"inq_zjNu_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSu",zjSu_ID)   ; call erreur(status,.TRUE.,"inq_zjSu_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEu",ziEu_ID)   ; call erreur(status,.TRUE.,"inq_ziEu_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWu",ziWu_ID)   ; call erreur(status,.TRUE.,"inq_ziWu_ID")
!-
status = NF90_INQ_VARID(fidcoeff,"wgNEv",wgNEv_ID) ; call erreur(status,.TRUE.,"inq_wgNEv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgNWv",wgNWv_ID) ; call erreur(status,.TRUE.,"inq_wgNWv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSWv",wgSWv_ID) ; call erreur(status,.TRUE.,"inq_wgSWv_ID")
status = NF90_INQ_VARID(fidcoeff,"wgSEv",wgSEv_ID) ; call erreur(status,.TRUE.,"inq_wgSEv_ID")
status = NF90_INQ_VARID(fidcoeff,"zjNv",zjNv_ID)   ; call erreur(status,.TRUE.,"inq_zjNv_ID")
status = NF90_INQ_VARID(fidcoeff,"zjSv",zjSv_ID)   ; call erreur(status,.TRUE.,"inq_zjSv_ID")
status = NF90_INQ_VARID(fidcoeff,"ziEv",ziEv_ID)   ; call erreur(status,.TRUE.,"inq_ziEv_ID")
status = NF90_INQ_VARID(fidcoeff,"ziWv",ziWv_ID)   ; call erreur(status,.TRUE.,"inq_ziWv_ID")
          
status = NF90_GET_VAR(fidcoeff,wgNEt_ID,wgNEt)     ; call erreur(status,.TRUE.,"getvar_wgNEt")
status = NF90_GET_VAR(fidcoeff,wgNWt_ID,wgNWt)     ; call erreur(status,.TRUE.,"getvar_wgNWt")
status = NF90_GET_VAR(fidcoeff,wgSWt_ID,wgSWt)     ; call erreur(status,.TRUE.,"getvar_wgSWt")
status = NF90_GET_VAR(fidcoeff,wgSEt_ID,wgSEt)     ; call erreur(status,.TRUE.,"getvar_wgSEt")
status = NF90_GET_VAR(fidcoeff,zjNt_ID,zjNt)       ; call erreur(status,.TRUE.,"getvar_zjNt")
status = NF90_GET_VAR(fidcoeff,zjSt_ID,zjSt)       ; call erreur(status,.TRUE.,"getvar_zjSt")
status = NF90_GET_VAR(fidcoeff,ziEt_ID,ziEt)       ; call erreur(status,.TRUE.,"getvar_ziEt")
status = NF90_GET_VAR(fidcoeff,ziWt_ID,ziWt)       ; call erreur(status,.TRUE.,"getvar_ziWt")
!-
status = NF90_GET_VAR(fidcoeff,wgNEu_ID,wgNEu)     ; call erreur(status,.TRUE.,"getvar_wgNEu")
status = NF90_GET_VAR(fidcoeff,wgNWu_ID,wgNWu)     ; call erreur(status,.TRUE.,"getvar_wgNWu")
status = NF90_GET_VAR(fidcoeff,wgSWu_ID,wgSWu)     ; call erreur(status,.TRUE.,"getvar_wgSWu")
status = NF90_GET_VAR(fidcoeff,wgSEu_ID,wgSEu)     ; call erreur(status,.TRUE.,"getvar_wgSEu")
status = NF90_GET_VAR(fidcoeff,zjNu_ID,zjNu)       ; call erreur(status,.TRUE.,"getvar_zjNu")
status = NF90_GET_VAR(fidcoeff,zjSu_ID,zjSu)       ; call erreur(status,.TRUE.,"getvar_zjSu")
status = NF90_GET_VAR(fidcoeff,ziEu_ID,ziEu)       ; call erreur(status,.TRUE.,"getvar_ziEu")
status = NF90_GET_VAR(fidcoeff,ziWu_ID,ziWu)       ; call erreur(status,.TRUE.,"getvar_ziWu")
!-
status = NF90_GET_VAR(fidcoeff,wgNEv_ID,wgNEv)     ; call erreur(status,.TRUE.,"getvar_wgNEv")
status = NF90_GET_VAR(fidcoeff,wgNWv_ID,wgNWv)     ; call erreur(status,.TRUE.,"getvar_wgNWv")
status = NF90_GET_VAR(fidcoeff,wgSWv_ID,wgSWv)     ; call erreur(status,.TRUE.,"getvar_wgSWv")
status = NF90_GET_VAR(fidcoeff,wgSEv_ID,wgSEv)     ; call erreur(status,.TRUE.,"getvar_wgSEv")
status = NF90_GET_VAR(fidcoeff,zjNv_ID,zjNv)       ; call erreur(status,.TRUE.,"getvar_zjNv")
status = NF90_GET_VAR(fidcoeff,zjSv_ID,zjSv)       ; call erreur(status,.TRUE.,"getvar_zjSv")
status = NF90_GET_VAR(fidcoeff,ziEv_ID,ziEv)       ; call erreur(status,.TRUE.,"getvar_ziEv")
status = NF90_GET_VAR(fidcoeff,ziWv_ID,ziWv)       ; call erreur(status,.TRUE.,"getvar_ziWv")
                            
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
          write(file_in_SSH,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_ssh)
        CASE(192)
          write(file_in_SSH,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_ssh) 
        CASE(193)
          write(file_in_SSH,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_ssh)
        CASE(194)
          write(file_in_SSH,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_ssh)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
     END SELECT
     write(*,*) 'Looking for existence of ', TRIM(file_in_SSH)
     inquire(file=file_in_SSH, exist=existfile)
     if ( existfile ) then; write(*,*) 'BINGO !'; exit ; endif
  enddo !-kfmt

  IF ( existfile ) THEN

    write(*,*) 'Reading T,S input dimensions in ', TRIM(file_in_SSH)
    status = NF90_OPEN(TRIM(file_in_SSH),0,fidT)            ; call erreur(status,.TRUE.,"read first")

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
    
    if ( mlon .ne. mx_GLO .or. mlat .ne. my_GLO .or. mdeptht .ne. mz_GLO ) then
      write(*,*) '~!@#$%^ MISMATCH BETWEEN DIMENSIONS OF GLOBAL FILES >>>>>>>>>>>> STOP !!'
      stop
    endif

    ALLOCATE( nav_lon(mx_GLO,my_GLO), nav_lat(mx_GLO,my_GLO) )
    ALLOCATE( deptht(mz_GLO) )

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
! 5- Process all SSH files over specified period
!=================================================================================

DO kyear=nn_yeari,nn_yearf

  DO kmonth=1,12

    DO kday=1,31

      SELECT CASE(nfmt)
        CASE(191)
          write(file_in_SSH,191) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_ssh)
        CASE(192)
          write(file_in_SSH,192) TRIM(data_dir), kyear, TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_ssh) 
        CASE(193)
          write(file_in_SSH,193) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(sep2), kday, TRIM(data_suffix_ssh)
        CASE(194)
          write(file_in_SSH,194) TRIM(data_dir), TRIM(data_prefix), kyear, TRIM(sep1), kmonth, TRIM(data_suffix_ssh)
        CASE DEFAULT 
          write(*,*) 'wrong nfmt value >>>>>> stop !'
          stop
      END SELECT
      inquire(file=file_in_SSH, exist=existfile)

      IF ( existfile ) THEN

        ! output file format :
        if     ( nfmt .eq. 191 .or. nfmt .eq. 193 ) then
          401 FORMAT(a,'/BDY/bdyT_ssh_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          411 FORMAT(a,'/BDY/bdyU_ssh_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          421 FORMAT(a,'/BDY/bdyV_ssh_',i4.4,'_',i2.2,'_',i2.2,'_',a,'.nc')
          write(file_bdy_SSH,  401) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          write(file_bdy_SSH_U,411) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
          write(file_bdy_SSH_V,421) TRIM(config_dir), kyear, kmonth, kday, TRIM(config)
        elseif ( nfmt .eq. 192 .or. nfmt .eq. 194 ) then
          402 FORMAT(a,'/BDY/bdyT_ssh_',i4.4,'_',i2.2,'_',a,'.nc')
          412 FORMAT(a,'/BDY/bdyU_ssh_',i4.4,'_',i2.2,'_',a,'.nc')
          422 FORMAT(a,'/BDY/bdyV_ssh_',i4.4,'_',i2.2,'_',a,'.nc')
          write(file_bdy_SSH,  402) TRIM(config_dir), kyear, kmonth, TRIM(config)
          write(file_bdy_SSH_U,412) TRIM(config_dir), kyear, kmonth, TRIM(config)
          write(file_bdy_SSH_V,422) TRIM(config_dir), kyear, kmonth, TRIM(config)
        else
          write(*,*) 'Do not forget to include new file format in the format definition for file_bdy_SSH  >>>> stop'
          stop
        endif

        ALLOCATE( ssh_T_GLO(mx_GLO,my_GLO,mtime)  )
        ALLOCATE( ssh_U_GLO(mx_GLO,my_GLO,mtime)  )
        ALLOCATE( ssh_V_GLO(mx_GLO,my_GLO,mtime)  )
        ALLOCATE( time(mtime) )
        
        !---------------------------------------
        ! 5a- Read input (GLOBAL) ssh :

        write(*,*) 'Reading ssh in ', TRIM(file_in_SSH)
        
        status = NF90_OPEN(TRIM(file_in_SSH),0,fidT)                  ; call erreur(status,.TRUE.,"read ORCA12 TS") 
        
        status = NF90_INQ_VARID(fidT,"time_counter",time_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"time",time_ID)
        call erreur(status,.TRUE.,"inq_time_ID")
        status = NF90_INQ_VARID(fidT,"sossheig",ssh_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"ssh",ssh_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"SSH",ssh_ID)
        if ( status .ne. 0 ) status = NF90_INQ_VARID(fidT,"zos",ssh_ID)
        call erreur(status,.TRUE.,"inq_ssh_ID")
        
        status = NF90_GET_VAR(fidT,time_ID,time)                      ; call erreur(status,.TRUE.,"getvar_time")
        status = NF90_GET_VAR(fidT,ssh_ID,ssh_T_GLO)                ; call erreur(status,.TRUE.,"getvar_ssh")

        status = NF90_GET_ATT(fidT,time_ID,"calendar",calendar)       ; call erreur(status,.TRUE.,"getatt_origin")
        status = NF90_GET_ATT(fidT,time_ID,"units",time_units)        ; call erreur(status,.TRUE.,"getatt_units")
        
        status = NF90_CLOSE(fidT)                                     ; call erreur(status,.TRUE.,"fin_lecture")     

        !---------------------------------------
        ! 5b- Calculate SSH on U and V grids :

        do l=1,mtime

          ssh_U_GLO(mx_GLO,:,l) = ssh_T_GLO(mx_GLO,:,l) ! to improve if periodic grid
          ssh_U_GLO(1:mx_GLO-1,:,l) = 0.5 * (   ssh_T_GLO(1:mx_GLO-1,:,l) * areaT(1:mx_GLO-1,:)    &
          &                                   + ssh_T_GLO(2:mx_GLO  ,:,l) * areaT(2:mx_GLO  ,:)  ) &
          &                                 / (eps + areaU(1:mx_GLO-1,:) )

          ssh_V_GLO(:,my_GLO,l) = ssh_T_GLO(:,my_GLO,l)
          ssh_V_GLO(:,1:my_GLO-1,l) = 0.5 * (   ssh_T_GLO(:,1:my_GLO-1,l) * areaT(:,1:my_GLO-1)    &
          &                                   + ssh_T_GLO(:,2:my_GLO  ,l) * areaT(:,2:my_GLO  )  ) &
          &                                 / (eps + areaV(:,1:my_GLO-1) )

        enddo

        !---------------------------------------
        ! 5c- Remove possible NaNs :
 
        do i=1,mx_GLO
        do j=1,my_GLO
        do l=1,mtime
          if ( .not. abs(ssh_T_GLO(i,j,l)) .lt. 100.0 ) ssh_T_GLO(i,j,l) = 0.0
          if ( .not. abs(ssh_U_GLO(i,j,l)) .lt. 100.0 ) ssh_U_GLO(i,j,l) = 0.0
          if ( .not. abs(ssh_V_GLO(i,j,l)) .lt. 100.0 ) ssh_V_GLO(i,j,l) = 0.0
        enddo
        enddo
        enddo

        !---------------------------------------
        ! 5d- Calculate SSH values on bdyT :
      
        ALLOCATE( ssh_T_bdy(mxbt,1,mtime)  )
 
        do kbdy=1,mxbt
           aNW = wgNWt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)) ,1 )
           aSW = wgSWt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)) ,1 )
           aNE = wgNEt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)) ,1 )
           aSE = wgSEt(nbit(kbdy,1),nbjt(kbdy,1)) * tmask_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)) ,1 )
           aa = aNW + aSW + aNE + aSE
           if ( aa .gt. eps .and. zjSt(nbit(kbdy,1),nbjt(kbdy,1)) .gt. 1 ) then
             do kt=1,mtime
                ssh_T_bdy(kbdy,1,kt) = (   ssh_T_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kt ) * aNW   &
                &                        + ssh_T_GLO( ziWt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kt ) * aSW   &
                &                        + ssh_T_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjNt(nbit(kbdy,1),nbjt(kbdy,1)), kt ) * aNE   &
                &                        + ssh_T_GLO( ziEt(nbit(kbdy,1),nbjt(kbdy,1)), zjSt(nbit(kbdy,1),nbjt(kbdy,1)), kt ) * aSE ) &
                &                           * tmask_REG(nbit(kbdy,1),nbjt(kbdy,1),1) / aa
             enddo
           else
             ssh_T_bdy(kbdy,1,:) = 0.e0
           endif
        enddo

        !---------------------------------------
        ! 5e- Calculate SSH values on bdyU :
      
        ALLOCATE( ssh_U_bdy(mxbu,1,mtime)  )
 
        do kbdy=1,mxbu
           aNW = wgNWu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)) ,1 )
           aSW = wgSWu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)) ,1 )
           aNE = wgNEu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)) ,1 )
           aSE = wgSEu(nbiu(kbdy,1),nbju(kbdy,1)) * umask_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)) ,1 )
           aa = aNW + aSW + aNE + aSE
           if ( aa .gt. eps .and. zjSu(nbiu(kbdy,1),nbju(kbdy,1)) .gt. 1 ) then
             do kt=1,mtime
                ssh_U_bdy(kbdy,1,kt) = (   ssh_U_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kt ) * aNW   &
                &                        + ssh_U_GLO( ziWu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kt ) * aSW   &
                &                        + ssh_U_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjNu(nbiu(kbdy,1),nbju(kbdy,1)), kt ) * aNE   &
                &                        + ssh_U_GLO( ziEu(nbiu(kbdy,1),nbju(kbdy,1)), zjSu(nbiu(kbdy,1),nbju(kbdy,1)), kt ) * aSE ) &
                &                           * umask_REG(nbiu(kbdy,1),nbju(kbdy,1),1) / aa
             enddo
           else
             ssh_U_bdy(kbdy,1,:) = 0.e0
           endif
        enddo

        !---------------------------------------
        ! 5f- Calculate SSH values on bdyV :
      
        ALLOCATE( ssh_V_bdy(mxbv,1,mtime)  )
 
        do kbdy=1,mxbv
           aNW = wgNWv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)) ,1 )
           aSW = wgSWv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)) ,1 )
           aNE = wgNEv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)) ,1 )
           aSE = wgSEv(nbiv(kbdy,1),nbjv(kbdy,1)) * vmask_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)) ,1 )
           aa = aNW + aSW + aNE + aSE
           if ( aa .gt. eps .and. zjSv(nbiv(kbdy,1),nbjv(kbdy,1)) .gt. 1 ) then
             do kt=1,mtime
                ssh_V_bdy(kbdy,1,kt) = (   ssh_V_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kt ) * aNW   &
                &                        + ssh_V_GLO( ziWv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kt ) * aSW   &
                &                        + ssh_V_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjNv(nbiv(kbdy,1),nbjv(kbdy,1)), kt ) * aNE   &
                &                        + ssh_V_GLO( ziEv(nbiv(kbdy,1),nbjv(kbdy,1)), zjSv(nbiv(kbdy,1),nbjv(kbdy,1)), kt ) * aSE ) &
                &                           * vmask_REG(nbiv(kbdy,1),nbjv(kbdy,1),1) / aa
             enddo
           else
             ssh_V_bdy(kbdy,1,:) = 0.e0
           endif
        enddo

        !---------------------------------------------------------
        ! 5g- Write new BDY netcdf file containing SSH at T points :

        write(*,*) 'Creating ', TRIM(file_bdy_SSH)
        status = NF90_CREATE(TRIM(file_bdy_SSH),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbT",mxbt,dimID_xbT)                             ; call erreur(status,.TRUE.,"def_dimID_xbT")

        status = NF90_DEF_VAR(fidM,"sossheig",NF90_FLOAT,(/dimID_xbT,dimID_yb,dimID_time_counter/),ssh_ID)
        call erreur(status,.TRUE.,"def_var_ssh_ID")
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
       
        status = NF90_PUT_ATT(fidM,ssh_ID,"long_name","Sea Surface Height")     ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,ssh_ID,"units","m")                          ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")          ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
        status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")           ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))    ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))   ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"long_name","bdy discrete distance") ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbrt_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbrt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"long_name","bdy j index")           ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbjt_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbjt_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"long_name","bdy i index")           ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        status = NF90_PUT_ATT(fidM,nbit_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbit_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_ssh.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,ssh_ID,ssh_T_bdy)     ; call erreur(status,.TRUE.,"var_ssh_ID")
        status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_bdy) ; call erreur(status,.TRUE.,"var_nav_lat_ID")
        status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_bdy) ; call erreur(status,.TRUE.,"var_nav_lon_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time) ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrt_ID,nbrt)         ; call erreur(status,.TRUE.,"var_nbrt_ID")
        status = NF90_PUT_VAR(fidM,nbjt_ID,nbjt)         ; call erreur(status,.TRUE.,"var_nbjt_ID")
        status = NF90_PUT_VAR(fidM,nbit_ID,nbit)         ; call erreur(status,.TRUE.,"var_nbit_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")

        !---------------------------------------------------------
        ! 5h- Write new BDY netcdf file containing SSH at U points :

        write(*,*) 'Creating ', TRIM(file_bdy_SSH_U)
        status = NF90_CREATE(TRIM(file_bdy_SSH_U),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbU",mxbu,dimID_xbU)                             ; call erreur(status,.TRUE.,"def_dimID_xbU")

        status = NF90_DEF_VAR(fidM,"ssh_U",NF90_FLOAT,(/dimID_xbU,dimID_yb,dimID_time_counter/),ssh_ID)
        call erreur(status,.TRUE.,"def_var_ssh_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbru_ID)
        call erreur(status,.TRUE.,"def_var_nbru_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbU,dimID_yb/),nbju_ID)
        call erreur(status,.TRUE.,"def_var_nbju_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbU,dimID_yb/),nbiu_ID)
        call erreur(status,.TRUE.,"def_var_nbiu_ID")
       
        status = NF90_PUT_ATT(fidM,ssh_ID,"long_name","SSH on U grid")          ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,ssh_ID,"units","m")                          ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))    ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))   ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbru_ID,"long_name","bdy discrete distance") ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidM,nbru_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbru_ID")
        status = NF90_PUT_ATT(fidM,nbju_ID,"long_name","bdy j index")           ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidM,nbju_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbju_ID")
        status = NF90_PUT_ATT(fidM,nbiu_ID,"long_name","bdy i index")           ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        status = NF90_PUT_ATT(fidM,nbiu_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbiu_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_ssh.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,ssh_ID,ssh_U_bdy)     ; call erreur(status,.TRUE.,"var_ssh_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time) ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbru_ID,nbru)         ; call erreur(status,.TRUE.,"var_nbru_ID")
        status = NF90_PUT_VAR(fidM,nbju_ID,nbju)         ; call erreur(status,.TRUE.,"var_nbju_ID")
        status = NF90_PUT_VAR(fidM,nbiu_ID,nbiu)         ; call erreur(status,.TRUE.,"var_nbiu_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")

        !---------------------------------------------------------
        ! 5i- Write new BDY netcdf file containing SSH at V points :

        write(*,*) 'Creating ', TRIM(file_bdy_SSH_V)
        status = NF90_CREATE(TRIM(file_bdy_SSH_V),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create BDY file')                     

        status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
        status = NF90_DEF_DIM(fidM,"yb",myb,dimID_yb)                                ; call erreur(status,.TRUE.,"def_dimID_yb")
        status = NF90_DEF_DIM(fidM,"xbV",mxbv,dimID_xbV)                             ; call erreur(status,.TRUE.,"def_dimID_xbV")

        status = NF90_DEF_VAR(fidM,"ssh_V",NF90_FLOAT,(/dimID_xbV,dimID_yb,dimID_time_counter/),ssh_ID)
        call erreur(status,.TRUE.,"def_var_ssh_ID")
        status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID)
        call erreur(status,.TRUE.,"def_var_time_counter_ID")
        status = NF90_DEF_VAR(fidM,"nbrdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbrv_ID)
        call erreur(status,.TRUE.,"def_var_nbrv_ID")
        status = NF90_DEF_VAR(fidM,"nbjdta",NF90_INT,(/dimID_xbV,dimID_yb/),nbjv_ID)
        call erreur(status,.TRUE.,"def_var_nbjv_ID")
        status = NF90_DEF_VAR(fidM,"nbidta",NF90_INT,(/dimID_xbV,dimID_yb/),nbiv_ID)
        call erreur(status,.TRUE.,"def_var_nbiv_ID")
       
        status = NF90_PUT_ATT(fidM,ssh_ID,"long_name","SSH on V grid")          ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,ssh_ID,"units","m")                          ; call erreur(status,.TRUE.,"put_att_ssh_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"units",TRIM(time_units))    ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar",TRIM(calendar))   ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
        status = NF90_PUT_ATT(fidM,nbrv_ID,"long_name","bdy discrete distance") ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidM,nbrv_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbrv_ID")
        status = NF90_PUT_ATT(fidM,nbjv_ID,"long_name","bdy j index")           ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidM,nbjv_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbjv_ID")
        status = NF90_PUT_ATT(fidM,nbiv_ID,"long_name","bdy i index")           ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        status = NF90_PUT_ATT(fidM,nbiv_ID,"units","unitless")                  ; call erreur(status,.TRUE.,"put_att_nbiv_ID")
        
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bdy_ssh.f90")
        status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
        call erreur(status,.TRUE.,"put_att_GLOBAL")
        
        status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"end definition") 
        
        status = NF90_PUT_VAR(fidM,ssh_ID,ssh_V_bdy)     ; call erreur(status,.TRUE.,"var_ssh_ID")
        status = NF90_PUT_VAR(fidM,time_counter_ID,time) ; call erreur(status,.TRUE.,"var_time_counter_ID")
        status = NF90_PUT_VAR(fidM,nbrv_ID,nbrv)         ; call erreur(status,.TRUE.,"var_nbrv_ID")
        status = NF90_PUT_VAR(fidM,nbjv_ID,nbjv)         ; call erreur(status,.TRUE.,"var_nbjv_ID")
        status = NF90_PUT_VAR(fidM,nbiv_ID,nbiv)         ; call erreur(status,.TRUE.,"var_nbiv_ID")
        
        status = NF90_CLOSE(fidM) ; call erreur(status,.TRUE.,"close BDY file")
        
        !---------------------------------------------------------

        DEALLOCATE( ssh_T_GLO, ssh_U_GLO, ssh_V_GLO, time )
        DEALLOCATE( ssh_T_bdy, ssh_U_bdy, ssh_V_bdy )

        !---------------------------------------------------------

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
