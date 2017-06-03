program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, Jan. 2017
!
! Script to extract the initial state from a global grid (e.g. eORCA12, ORCA025).
!
! 1- Read GLOBAL mask
! 2- Read REGIONAL mask
! 3- Read GLOBAL temperature initial state
! 4- Read GLOBAL salinity initial state
! 5- Extract regional from global
! 6- Writing initial state for temperature
! 7- Writing initial state for salinity
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
namelist /init/ nn_init, file_in_mask_extract, file_in_T, file_in_S, nn_eosmatch, nn_iter, nn_rsmax, nn_rzmax, &
&               rn_temp, rn_sal, nn_smooth, file_in_zgr_extract
INTEGER                               :: nn_init, nn_iter, nn_rsmax, nn_rzmax, nn_eosmatch, nn_smooth
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: file_in_mask_extract, file_in_zgr_extract, config_dir, file_in_T, file_in_S
REAL(KIND=4)                          :: rn_temp, rn_sal

INTEGER :: fidMSKIN, fidMSKREG, status, dimID_z, dimID_y, dimID_x, mz_GLO, my_GLO, mx_GLO, tmask_GLO_ID, &
&          mx_REG, my_REG, mz_REG, tmask_REG_ID, fidSAL, fidTEMP, votemper_ID, vosaline_ID,              &
&          dimID_time_counter, ai, bi, aj, bj, iii, jjj, kkk, kk, iGLO, jGLO, iREG, jREG, fidTin, fidSin,&
&          kiter, rs, rz, sg, time_counter_ID, fidCOORD, imin_ORCA12, jmin_ORCA12, lon_ID, lat_ID, dij,  &
&          dep_ID, kGLO, ntest, im1, ip1, jm1, jp1, mx_tmp, my_tmp, fidZGR, e3t_GLO_ID

CHARACTER(LEN=180) :: file_in_mask_REG, file_in_coord_REG, file_out_temp, file_out_sal 

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: tmask_GLO, tmask_REG, missing, tmp_missing

REAL(KIND=4),ALLOCATABLE,DIMENSION(:) ::  dep_GLO

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: lon_GLO, lat_GLO

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:) :: votemper_GLO, vosaline_GLO, votemper_REG, vosaline_REG, &
&                                            tmp_votemper_REG, tmp_vosaline_REG, e3t_GLO

LOGICAL :: iout

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEt_ID, wgNWt_ID, wgSWt_ID, wgSEt_ID, angYt_ID, angXt_ID,  &
&                                         zkUt_ID, ziWt_ID, ziEt_ID, zjNt_ID, zjSt_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWt, ziEt, zjNt, zjSt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYt, angXt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEt, wgNWt, wgSWt, wgSEt

!-- Regional initial state
REAL*8                                 :: eps, aa, aSW, aNW, aSE, aNE

!=================================================================================
! 0- Initializations 
!=================================================================================

call gsw_saar_init (.true.)

! Default values (replaced with namelist values if specified):
config_dir        = '.'
nn_iter           = 100
nn_rsmax          =   5
nn_rzmax          =   1
nn_eosmatch       =   1
nn_smooth         =   1

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=init)
CLOSE(1)

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!- output file names :
write(file_out_temp,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/dta_temp_',a,'.nc')
write(file_out_sal,202)  TRIM(config_dir), TRIM(config)
202 FORMAT(a,'/dta_sal_',a,'.nc')

eps=1.d-9

! name of interpolation weights file :
write(file_coeff,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coeff_linear_',a,'.nc')

!=================================================================================
! 1a- Read GLOBAL mask :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_in_mask_extract),0,fidMSKIN); call erreur(status,.TRUE.,"read mask input") 

status = NF90_INQ_DIMID(fidMSKIN,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z_GLO")
status = NF90_INQ_DIMID(fidMSKIN,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_GLO")
status = NF90_INQ_DIMID(fidMSKIN,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_GLO")

status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_z,len=mz_GLO); call erreur(status,.TRUE.,"inq_dim_z_GLO")
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_y,len=my_GLO); call erreur(status,.TRUE.,"inq_dim_y_GLO")
status = NF90_INQUIRE_DIMENSION(fidMSKIN,dimID_x,len=mx_GLO); call erreur(status,.TRUE.,"inq_dim_x_GLO")

ALLOCATE(  tmask_GLO(mx_GLO,my_GLO,mz_GLO)  ) 
ALLOCATE(  lon_GLO  (mx_GLO,my_GLO)  )
ALLOCATE(  lat_GLO  (mx_GLO,my_GLO)  )
ALLOCATE(  dep_GLO  (mz_GLO)  )

status = NF90_INQ_VARID(fidMSKIN,"tmask",tmask_GLO_ID); call erreur(status,.TRUE.,"inq_tmask_GLO_ID")
status = NF90_INQ_VARID(fidMSKIN,"nav_lon",lon_ID)    ; call erreur(status,.TRUE.,"inq_lon_GLO_ID")
status = NF90_INQ_VARID(fidMSKIN,"nav_lat",lat_ID)    ; call erreur(status,.TRUE.,"inq_lat_GLO_ID")
status = NF90_INQ_VARID(fidMSKIN,"nav_lev",dep_ID)    ; call erreur(status,.TRUE.,"inq_dep_GLO_ID")

status = NF90_GET_VAR(fidMSKIN,tmask_GLO_ID,tmask_GLO); call erreur(status,.TRUE.,"getvar_tmask_GLO")
status = NF90_GET_VAR(fidMSKIN,lon_ID,lon_GLO)        ; call erreur(status,.TRUE.,"getvar_lon_GLO")
status = NF90_GET_VAR(fidMSKIN,lat_ID,lat_GLO)        ; call erreur(status,.TRUE.,"getvar_lat_GLO")
status = NF90_GET_VAR(fidMSKIN,dep_ID,dep_GLO)        ; call erreur(status,.TRUE.,"getvar_dep_GLO")

status = NF90_CLOSE(fidMSKIN); call erreur(status,.TRUE.,"end read mask_GLO")

!=================================================================================
! 1b- Read GLOBAL zgr (e3t) :                                 
!=================================================================================

status = NF90_OPEN(TRIM(file_in_zgr_extract),0,fidZGR); call erreur(status,.TRUE.,"read e3t_GLO") 

ALLOCATE(  e3t_GLO(mx_GLO,my_GLO,mz_GLO)  ) 

status = NF90_INQ_VARID(fidZGR,"e3t",e3t_GLO_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidZGR,"e3t_0",e3t_GLO_ID)
call erreur(status,.TRUE.,"inq_e3t_GLO_ID")

status = NF90_GET_VAR(fidZGR,e3t_GLO_ID,e3t_GLO); call erreur(status,.TRUE.,"getvar_e3t_GLO")

status = NF90_CLOSE(fidZGR); call erreur(status,.TRUE.,"end read e3t_GLO")

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
                                                     
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")

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
! 4- Read GLOBAL temperature initial state
!=================================================================================

status = NF90_OPEN(TRIM(file_in_T),0,fidTin); call erreur(status,.TRUE.,"read_T")        
ALLOCATE(  votemper_GLO(mx_GLO,my_GLO,mz_GLO)  )   
status = NF90_INQ_VARID(fidTin,"votemper",votemper_ID); call erreur(status,.TRUE.,"inq_votemper_GLO_ID")
status = NF90_GET_VAR(fidTin,votemper_ID,votemper_GLO); call erreur(status,.TRUE.,"getvar_votemper_GLO")
status = NF90_CLOSE(fidTin); call erreur(status,.TRUE.,"end_read_T") 

!=================================================================================
! 4- Read GLOBAL salinity initial state
!=================================================================================

status = NF90_OPEN(TRIM(file_in_S),0,fidSin); call erreur(status,.TRUE.,"read_S")
ALLOCATE(  vosaline_GLO(mx_GLO,my_GLO,mz_GLO)  )
status = NF90_INQ_VARID(fidSin,"vosaline",vosaline_ID); call erreur(status,.TRUE.,"inq_vosaline_GLO_ID")
status = NF90_GET_VAR(fidSin,vosaline_ID,vosaline_GLO); call erreur(status,.TRUE.,"getvar_vosaline_GLO")
status = NF90_CLOSE(fidSin); call erreur(status,.TRUE.,"end_read_S") 

!- convert to conservative temperature if needed :
if ( nn_eosmatch .eq. 0 ) then
  write(*,*) 'Converting from EOS80 to TEOS10 ...'
  do iGLO=1,mx_GLO
  do jGLO=1,my_GLO
  do kGLO=1,mz_GLO
    if ( tmask_GLO(iGLO,jGLO,kGLO) .eq. 1 ) then
      vosaline_GLO(iGLO,jGLO,kGLO) = gsw_sa_from_sp( DBLE(vosaline_GLO(iGLO,jGLO,kGLO)), DBLE(dep_GLO(kGLO)), DBLE(lon_GLO(iGLO,jGLO)), DBLE(lat_GLO(iGLO,jGLO)) )
      votemper_GLO(iGLO,jGLO,kGLO) = gsw_ct_from_pt( DBLE(vosaline_GLO(iGLO,jGLO,kGLO)), DBLE(votemper_GLO(iGLO,jGLO,kGLO)) )
    else
      vosaline_GLO(iGLO,jGLO,kGLO) = 0.d0
      votemper_GLO(iGLO,jGLO,kGLO) = 0.d0
    endif
    if ( votemper_GLO(iGLO,jGLO,kGLO) .lt. -2.5 ) then
      write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) iGLO,jGLO,kGLO
      write(*,*) lon_GLO(iGLO,jGLO), lat_GLO(iGLO,jGLO), dep_GLO(kGLO)
      write(*,*) vosaline_GLO(iGLO,jGLO,kGLO), votemper_GLO(iGLO,jGLO,kGLO)
    endif
  enddo
  enddo
  enddo
elseif ( nn_eosmatch .ne. 1 ) then
  write(*,*) '~!@#$%^* Error: nn_eosmatch should be 0 or 1 >>>>> stop !!'
  stop
endif

!=================================================================================
! 5- Extract regional from global
!=================================================================================

ALLOCATE( votemper_REG(mx_REG,my_REG,mz_REG) )
ALLOCATE( vosaline_REG(mx_REG,my_REG,mz_REG) )

IF ( nn_init == 1 ) THEN

  write(*,*) ' nn_init=1 HAS NOT BEEN IMPLEMENTED YET  >>>>> STOP'
  stop

ELSEIF ( nn_init == 2 ) THEN

  if ( mz_REG .ne. mz_GLO ) then
    write(*,*) '~!@#$%^ Adapt script for different number of vertical levels >>>> Stop!!'
    stop
  endif

 ! !- Read global attributes of coordinate file to get grid correspondance :
 ! !       i_ORCA12 = ai * i_ORCA025 + bi
 ! !       j_ORCA12 = aj * j_ORCA025 + bj
 ! status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidCOORD); call erreur(status,.TRUE.,"read coord input")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "ai", ai); call erreur(status,.TRUE.,"read att1")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "bi", bi); call erreur(status,.TRUE.,"read att2")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "aj", aj); call erreur(status,.TRUE.,"read att3")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "bj", bj); call erreur(status,.TRUE.,"read att4")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "imin_extraction", imin_ORCA12); call erreur(status,.TRUE.,"read att5")
 ! status = NF90_GET_ATT(fidCOORD, NF90_GLOBAL, "jmin_extraction", jmin_ORCA12); call erreur(status,.TRUE.,"read att6")
 ! status = NF90_CLOSE(fidCOORD)                         ; call erreur(status,.TRUE.,"end read fidCOORD")

  ! Just extract where ocean points on both grids :
  ALLOCATE( missing(mx_REG,my_REG,mz_REG) )
  ALLOCATE( tmp_missing(mx_REG,my_REG,mz_REG) )
  ALLOCATE( tmp_votemper_REG(mx_REG,my_REG,mz_REG) )
  ALLOCATE( tmp_vosaline_REG(mx_REG,my_REG,mz_REG) )
  missing(:,:,:)=0
  votemper_REG(:,:,:)=0.d0
  vosaline_REG(:,:,:)=0.d0
  do iREG=1,mx_REG
  do jREG=1,my_REG
  do kk=1,mz_REG

      aNW = tmask_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG),kk ) * e3t_GLO( ziWt(iREG,jREG),zjNt(iREG,jREG),kk ) * wgNWt(iREG,jREG)
      aSW = tmask_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG),kk ) * e3t_GLO( ziWt(iREG,jREG),zjSt(iREG,jREG),kk ) * wgSWt(iREG,jREG)
      aNE = tmask_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG),kk ) * e3t_GLO( ziEt(iREG,jREG),zjNt(iREG,jREG),kk ) * wgNEt(iREG,jREG)
      aSE = tmask_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG),kk ) * e3t_GLO( ziEt(iREG,jREG),zjSt(iREG,jREG),kk ) * wgSEt(iREG,jREG)

      aa = aNW + aSW + aNE + aSE

      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then

        votemper_REG (iREG,jREG,kk) = (   votemper_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG),kk ) * aNW       &
        &                               + votemper_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG),kk ) * aSW       &
        &                               + votemper_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG),kk ) * aNE       &
        &                               + votemper_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG),kk ) * aSE ) / aa

        vosaline_REG (iREG,jREG,kk) = (   vosaline_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG),kk ) * aNW       &
        &                               + vosaline_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG),kk ) * aSW       &
        &                               + vosaline_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG),kk ) * aNE       &
        &                               + vosaline_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG),kk ) * aSE ) / aa

      elseif ( tmask_REG(iREG,jREG,kk) .eq. 1 ) then !- oceanic point on regional grid but all land points on global grid

        missing(iREG,jREG,kk) = 1

      endif

  enddo
  enddo
  enddo

  ! Look for closest neighbours where we have missing values:
  do kiter=1,nn_iter
    ntest = NINT(sum(sum(sum(FLOAT(missing),3),2),1))
    write(*,*) '  kiter = ', kiter
    write(*,*) '     nb of pts with missing value: ', ntest
    if ( ntest .eq. 0 ) exit
    tmp_votemper_REG(:,:,:)=votemper_REG(:,:,:)
    tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
    tmp_missing(:,:,:)=missing(:,:,:)
    do iREG=1,mx_REG
    do jREG=1,my_REG
    do kk=1,mz_REG
      if ( missing(iREG,jREG,kk) .eq. 1 ) then
        iout=.FALSE.
        do rz=0,nn_rzmax,1
        do sg=-1,1,2 ! to look above first, then below
          do rs=1,nn_rsmax,1
            iii=iREG               ; jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG) ! to look right above/below
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=iREG               ; jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=iREG               ; jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG) 
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG) 
            if ( tmask_REG(iii,jjj,kkk) .eq. 1 .and. missing(iii,jjj,kkk) .eq. 0 ) then
              iout=.TRUE.
              exit
            endif
          enddo !- rs
          if (iout) exit
        enddo !-sg
        if (iout) then
          tmp_missing(iREG,jREG,kk) = 0
          tmp_votemper_REG(iREG,jREG,kk) = votemper_REG(iii,jjj,kkk)
          tmp_vosaline_REG(iREG,jREG,kk) = vosaline_REG(iii,jjj,kkk)
          exit
        elseif ( rz .eq. nn_rzmax .and. kiter .eq. nn_iter ) then
          write(*,953) iREG, jREG, kk
          953 FORMAT(' >>> WARNING for point (',3I5,') --> filled with rn_temp and rn_sal (to avoid this, increase nn_rsmax and/or nn_rzmax and/or nn_iter)')
          tmp_missing(iREG,jREG,kk) = 0
          tmp_votemper_REG(iREG,jREG,kk) = rn_temp
          tmp_vosaline_REG(iREG,jREG,kk) = rn_sal
          exit
        endif
        enddo !-rz
      endif !-if ( missing(iREG,jREG,kk) .eq. 1 )
    enddo !- kk
    enddo !- jREG
    enddo !- iREG
    missing(:,:,:)=tmp_missing(:,:,:)
    votemper_REG(:,:,:)=tmp_votemper_REG(:,:,:)
    vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
  enddo !- kiter

  !- Smoothing :
  if ( nn_smooth .gt. 1 ) then
    write(*,*) 'Smoothing window width = ', nn_smooth
    dij=INT(nn_smooth*0.5)
    tmp_votemper_REG(:,:,:)=votemper_REG(:,:,:)
    tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
    do iREG=1,mx_REG
    do jREG=1,my_REG
    do kk=1,mz_REG
      im1=MAX(iREG-dij,1) ; ip1=MIN(iREG+dij,mx_REG) 
      jm1=MAX(jREG-dij,1) ; jp1=MIN(jREG+dij,my_REG)
      if ( tmask_REG(iREG,jREG,kk) .eq. 1 ) then 
        tmp_votemper_REG(iREG,jREG,kk) =   SUM( SUM( votemper_REG(im1:ip1,jm1:jp1,kk) * tmask_REG(im1:ip1,jm1:jp1,kk), 2), 1) &
        &                                / SUM( SUM(                             1.0  * tmask_REG(im1:ip1,jm1:jp1,kk), 2), 1)
        tmp_vosaline_REG(iREG,jREG,kk) =   SUM( SUM( vosaline_REG(im1:ip1,jm1:jp1,kk) * tmask_REG(im1:ip1,jm1:jp1,kk), 2), 1) &
        &                                / SUM( SUM(                             1.0  * tmask_REG(im1:ip1,jm1:jp1,kk), 2), 1)
      else
        tmp_votemper_REG(iREG,jREG,kk) = 0.d0
        tmp_vosaline_REG(iREG,jREG,kk) = 0.d0
      endif
    enddo
    enddo
    enddo
    votemper_REG(:,:,:)=tmp_votemper_REG(:,:,:)
    vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
  else
    write(*,*) 'No Smoothing'
  endif

  !- "Drowning", i.e. put closest value everywhere on the mask file to avoid issue if namdom is slightly changed :
  !  We just repeat the previous methodology, but for masked points
  write(*,*) 'Drowning, i.e. fill all masked points with closest neighbour'
  missing(:,:,:)=NINT(1-FLOAT(tmask_REG(:,:,:)))
  ! Look for closest neighbours where we have missing values:
  do kiter=1,nn_iter
    ntest = NINT(sum(sum(sum(FLOAT(missing),3),2),1))
    write(*,*) '  kiter = ', kiter
    write(*,*) '     remaining nb of masked points to fill: ', ntest
    if ( ntest .eq. 0 ) exit
    tmp_votemper_REG(:,:,:)=votemper_REG(:,:,:)
    tmp_vosaline_REG(:,:,:)=vosaline_REG(:,:,:)
    tmp_missing(:,:,:)=missing(:,:,:)
    do iREG=1,mx_REG
    do jREG=1,my_REG
    do kk=1,mz_REG
      if ( missing(iREG,jREG,kk) .eq. 1 ) then
        iout=.FALSE.
        do rz=0,nn_rzmax,1
        do sg=-1,1,2 ! to look above first, then below
          do rs=1,nn_rsmax,1
            iii=iREG               ; jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG) ! to look right above/below
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=jREG               ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=iREG               ; jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=iREG               ; jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MIN(iREG+rs,mx_REG); jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG)
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=MIN(jREG+rs,my_REG); kkk= MIN(MAX(kk+rz*sg,1),mz_REG) 
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
            iii=MAX(iREG-rs,1)     ; jjj=MAX(jREG-rs,1)     ; kkk= MIN(MAX(kk+rz*sg,1),mz_REG) 
            if ( missing(iii,jjj,kkk) .eq. 0 ) then ; iout=.TRUE. ; exit ; endif
          enddo !- rs
          if (iout) exit
        enddo !-sg
        if (iout) then
          tmp_missing(iREG,jREG,kk) = 0
          tmp_votemper_REG(iREG,jREG,kk) = votemper_REG(iii,jjj,kkk)
          tmp_vosaline_REG(iREG,jREG,kk) = vosaline_REG(iii,jjj,kkk)
          exit
        elseif ( rz .eq. nn_rzmax .and. kiter .eq. nn_iter ) then
          tmp_missing(iREG,jREG,kk) = 0
          tmp_votemper_REG(iREG,jREG,kk) = rn_temp
          tmp_vosaline_REG(iREG,jREG,kk) = rn_sal
          exit
        endif
        enddo !-rz
      endif !-if ( missing(iREG,jREG,kk) .eq. 1 )
    enddo !- kk
    enddo !- jREG
    enddo !- iREG
    missing(:,:,:)=tmp_missing(:,:,:)
    votemper_REG(:,:,:)=tmp_votemper_REG(:,:,:)
    vosaline_REG(:,:,:)=tmp_vosaline_REG(:,:,:)
  enddo !- kiter

  !--  
  DEALLOCATE( tmp_votemper_REG, tmp_vosaline_REG, missing )

ELSE

  write(*,*) ' THIS nn_init VALUE DOES NOT CORRESPOND TO SOMETHING KNOWN  >>>>> STOP'
  stop

ENDIF

!=================================================================================
! 6- Writing initial state for temperature 
!=================================================================================

write(*,*) 'Writing ', TRIM(file_out_temp)

status = NF90_CREATE(TRIM(file_out_temp),NF90_NOCLOBBER,fidTEMP) ; call erreur(status,.TRUE.,'create output temp')

status = NF90_DEF_DIM(fidTEMP,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidTEMP,"x",mx_REG,dimID_x)                               ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidTEMP,"y",my_REG,dimID_y)                               ; call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidTEMP,"z",mz_REG,dimID_z)                               ; call erreur(status,.TRUE.,"def_dimID_deptht")

status = NF90_DEF_VAR(fidTEMP,"time_counter",NF90_DOUBLE,(/dimID_time_counter/),time_counter_ID)
call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidTEMP,"votemper",NF90_FLOAT,(/dimID_x,dimID_y,dimID_z,dimID_time_counter/),votemper_ID)
call erreur(status,.TRUE.,"def_var_votemper_ID")

status = NF90_PUT_ATT(fidTEMP,votemper_ID,"associate","time_counter, z, y, x") ; call erreur(status,.TRUE.,"put_att_votemper_ID")
status = NF90_PUT_ATT(fidTEMP,votemper_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_votemper_ID")
status = NF90_PUT_ATT(fidTEMP,votemper_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_votemper_ID")
status = NF90_PUT_ATT(fidTEMP,votemper_ID,"units","degC")                      ; call erreur(status,.TRUE.,"put_att_votemper_ID")
!if ( nn_eosmatch .eq. 0 ) then
  status = NF90_PUT_ATT(fidTEMP,votemper_ID,"long_name","conservative temperature") ; call erreur(status,.TRUE.,"put_att_votemper_ID")
!else
!  status = NF90_PUT_ATT(fidTEMP,votemper_ID,"long_name","potential temperature")    ; call erreur(status,.TRUE.,"put_att_votemper_ID") 
!endif
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"title","Time")                  ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"long_name","Time axis")         ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"standard_name","time")          ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidTEMP,time_counter_ID,"axis","T")                      ; call erreur(status,.TRUE.,"put_att_time_counter_ID")

status = NF90_PUT_ATT(fidTEMP,NF90_GLOBAL,"history","Created using extract_istate.f90")
status = NF90_PUT_ATT(fidTEMP,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidTEMP) ; call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidTEMP,time_counter_ID,1.0)       ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidTEMP,votemper_ID,votemper_REG)  ; call erreur(status,.TRUE.,"var_votemper_ID")

status = NF90_CLOSE(fidTEMP) ; call erreur(status,.TRUE.,"final")         

!=================================================================================
! 7- Writing initial state for salinity
!=================================================================================

write(*,*) 'Writing ', TRIM(file_out_sal)

status = NF90_CREATE(TRIM(file_out_sal),NF90_NOCLOBBER,fidSAL) ; call erreur(status,.TRUE.,'create output sal')

status = NF90_DEF_DIM(fidSAL,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidSAL,"x",mx_REG,dimID_x)                               ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidSAL,"y",my_REG,dimID_y)                               ; call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidSAL,"z",mz_REG,dimID_z)                               ; call erreur(status,.TRUE.,"def_dimID_deptht")

status = NF90_DEF_VAR(fidSAL,"time_counter",NF90_DOUBLE,(/dimID_time_counter/),time_counter_ID)
call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidSAL,"vosaline",NF90_FLOAT,(/dimID_x,dimID_y,dimID_z,dimID_time_counter/),vosaline_ID)
call erreur(status,.TRUE.,"def_var_vosaline_ID")

status = NF90_PUT_ATT(fidSAL,vosaline_ID,"associate","time_counter, z, y, x") ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
status = NF90_PUT_ATT(fidSAL,vosaline_ID,"missing_value",0.)                  ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
status = NF90_PUT_ATT(fidSAL,vosaline_ID,"_FillValue",0.)                     ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
!if ( nn_eosmatch .eq. 0 ) then
!  status = NF90_PUT_ATT(fidSAL,vosaline_ID,"units","psu")                     ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
!  status = NF90_PUT_ATT(fidSAL,vosaline_ID,"long_name","practical salinity")  ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
!else
  status = NF90_PUT_ATT(fidSAL,vosaline_ID,"units","g/kg")                    ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
  status = NF90_PUT_ATT(fidSAL,vosaline_ID,"long_name","absolute salinity")   ; call erreur(status,.TRUE.,"put_att_vosaline_ID")
!endif
status = NF90_PUT_ATT(fidSAL,time_counter_ID,"title","Time")                  ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidSAL,time_counter_ID,"long_name","Time axis")         ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidSAL,time_counter_ID,"standard_name","time")          ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidSAL,time_counter_ID,"axis","T")                      ; call erreur(status,.TRUE.,"put_att_time_counter_ID")

status = NF90_PUT_ATT(fidSAL,NF90_GLOBAL,"history","Created using extract_istate.f90") ; call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidSAL) ; call erreur(status,.TRUE.,"fin_definition") 

status = NF90_PUT_VAR(fidSAL,time_counter_ID,1.0)       ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidSAL,vosaline_ID,vosaline_REG)  ; call erreur(status,.TRUE.,"var_vosaline_ID")

status = NF90_CLOSE(fidSAL) ; call erreur(status,.TRUE.,"final")         


end program modif



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
