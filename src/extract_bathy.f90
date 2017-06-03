program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, Jan. 2017
!
! Script to extract the bathymetry from a global grid (e.g. eORCA12).
!
! The bathymetry along the boundaries (over a 3-pts halo) is the same as in the
! dataset used as lateral boundary conditions.
!
! 0- Initializations 
! 1- Read coarse bathymetry used for consistent bathymetry along boundaries
! 2- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
! 3- Extract variables on the REG grid :
! 4- Manual corrections for WED12 :
! 5- Writing new REG bathymetry file :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_in_bathy_bdy, ln_isfcav, &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, file_in_coord_bdy
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir, file_in_coord_bdy
LOGICAL                               :: ln_isfcav

!-- local variables :
INTEGER :: fidORCA12, fidBDY, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_ORCA12, mx_ORCA12,  my_BDY, mx_BDY,  my_REG, mx_REG, imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12,            &
&          ai, aj, bi, bj, iREG, jREG, jtmp, npts, kk, ki, kj, ni1, ni2, nj1, nj2, pi, pj, kiref, kjref, mx_tmp, my_tmp, e2f_ID, e2v_ID,  &
&          e2u_ID, e2t_ID, e1f_ID, e1v_ID, e1u_ID, e1t_ID, gphif_ID, gphiv_ID, gphiu_ID, gphit_ID, glamf_ID, glamv_ID, glamu_ID, glamt_ID,&
&          fidCOORDreg, fidCOORDpar 

CHARACTER(LEN=150) :: aaa, file_bathy_out, file_in_coord_REG

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: nav_lat_GLO, nav_lon_GLO, isf_draft_GLO, Bathymetry_isf_GLO, Bathymetry_GLO,  &
&                                          nav_lat_BDY, nav_lon_BDY, Bathymetry_BDY,                                        &
&                                          nav_lat_REG, nav_lon_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: e2f, e2v, e2u, e2t, e1f, e1v, e1u, e1t, gphif, gphiv, gphiu, gphit, glamf, glamv, glamu, glamt, &
&                                          e2f_REG, e2v_REG, e2u_REG, e2t_REG, e1f_REG, e1v_REG, e1u_REG, e1t_REG,                         &
&                                          gphif_REG, gphiv_REG, gphiu_REG, gphit_REG, glamf_REG, glamv_REG, glamu_REG, glamt_REG

INTEGER(KIND=4), DIMENSION(8) :: imin, imax, jmin

!-- interpolation parameters
INTEGER                                :: fidcoeff, wgNEt_ID, wgNWt_ID, wgSWt_ID, wgSEt_ID, angYt_ID, angXt_ID,  &
&                                         zkUt_ID, ziWt_ID, ziEt_ID, zjNt_ID, zjSt_ID
CHARACTER(LEN=150)                     :: file_coeff
INTEGER*2,ALLOCATABLE,DIMENSION(:,:)   :: ziWt, ziEt, zjNt, zjSt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: angYt, angXt
REAL*4,ALLOCATABLE,DIMENSION(:,:)      :: wgNEt, wgNWt, wgSWt, wgSEt

!-- Regional initial state
REAL*8                                 :: eps, aa
 
!=================================================================================
! 0- Initializations 
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'
file_in_bathy_bdy = 'not_used'
ln_isfcav         = .false.

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
CLOSE(1)

! name of regional bathymetry file (output file) :
write(file_bathy_out,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/bathy_meter_',a,'.nc')

! name of regional coordinates file (output file) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

! name of interpolation weights file :
write(file_coeff,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coeff_linear_',a,'.nc')

!-
imin_ORCA12 = nn_imin_extract
imax_ORCA12 = nn_imax_extract
jmin_ORCA12 = nn_jmin_extract
jmax_ORCA12 = nn_jmax_extract

eps = 1.d-9

!=================================================================================
! 1a- Read grid correspondance with GLO (i.e. extraction coordinates)
!=================================================================================

status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidCOORDreg); call erreur(status,.TRUE.,"read coord input")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "imin_extraction", imin_ORCA12); call erreur(status,.TRUE.,"read att5")
status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "jmin_extraction", jmin_ORCA12); call erreur(status,.TRUE.,"read att6")
status = NF90_CLOSE(fidCOORDreg)                         ; call erreur(status,.TRUE.,"end read fidCOORDreg")

!=================================================================================
! 1b- Read bathymetry in which to extract (e.g. eORCA12)
!=================================================================================

write(*,*) 'Extracting regional bathymetry from : ', TRIM(file_in_bathy_extract)

status = NF90_OPEN(TRIM(file_in_bathy_extract),0,fidORCA12) ; call erreur(status,.TRUE.,"read_bathy_to_extract") 

status = NF90_INQ_DIMID(fidORCA12,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_ORCA12")
status = NF90_INQ_DIMID(fidORCA12,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_ORCA12")

status = NF90_INQUIRE_DIMENSION(fidORCA12,dimID_y,len=my_ORCA12); call erreur(status,.TRUE.,"inq_dim_y_ORCA12")
status = NF90_INQUIRE_DIMENSION(fidORCA12,dimID_x,len=mx_ORCA12); call erreur(status,.TRUE.,"inq_dim_x_ORCA12")

ALLOCATE(  nav_lat_GLO        (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  nav_lon_GLO        (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  isf_draft_GLO      (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  Bathymetry_isf_GLO (mx_ORCA12,my_ORCA12)  ) 
ALLOCATE(  Bathymetry_GLO     (mx_ORCA12,my_ORCA12)  ) 

status = NF90_INQ_VARID(fidORCA12,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_ORCA12")
status = NF90_INQ_VARID(fidORCA12,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidORCA12,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_ORCA12")
status = NF90_INQ_VARID(fidORCA12,"Bathymetry",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_ORCA12")

status = NF90_GET_VAR(fidORCA12,nav_lat_ID,nav_lat_GLO); call erreur(status,.TRUE.,"getvar_nav_lat_GLO")
status = NF90_GET_VAR(fidORCA12,nav_lon_ID,nav_lon_GLO); call erreur(status,.TRUE.,"getvar_nav_lon_GLO")
status = NF90_GET_VAR(fidORCA12,Bathymetry_ID,Bathymetry_GLO); call erreur(status,.TRUE.,"getvar_Bathymetry_GLO")

if ( ln_isfcav ) then
  status = NF90_INQ_VARID(fidORCA12,"isf_draft",isf_draft_ID);              call erreur(status,.TRUE.,"inq_isf_draft_ID_ORCA12")
  status = NF90_INQ_VARID(fidORCA12,"Bathymetry_isf",Bathymetry_isf_ID);    call erreur(status,.TRUE.,"inq_Bathymetry_isf_ID_ORCA12")
  status = NF90_GET_VAR(fidORCA12,isf_draft_ID,isf_draft_GLO);           call erreur(status,.TRUE.,"getvar_isf_draft_GLO")
  status = NF90_GET_VAR(fidORCA12,Bathymetry_isf_ID,Bathymetry_isf_GLO); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_GLO")
endif

status = NF90_CLOSE(fidORCA12); call erreur(status,.TRUE.,"close_grid_to_extract")     


!=================================================================================
! 1- Read coarse bathymetry used for consistent bathymetry along boundaries
!=================================================================================

write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)

status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidBDY); call erreur(status,.TRUE.,"read_coarse_bathymetry") 

status = NF90_INQ_DIMID(fidBDY,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_BDY")
status = NF90_INQ_DIMID(fidBDY,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_BDY")

status = NF90_INQUIRE_DIMENSION(fidBDY,dimID_y,len=my_BDY); call erreur(status,.TRUE.,"inq_dim_y_BDY")
status = NF90_INQUIRE_DIMENSION(fidBDY,dimID_x,len=mx_BDY); call erreur(status,.TRUE.,"inq_dim_x_BDY")

ALLOCATE(  Bathymetry_BDY (mx_BDY,my_BDY)  ) 
ALLOCATE(  nav_lat_BDY    (mx_BDY,my_BDY)  ) 
ALLOCATE(  nav_lon_BDY    (mx_BDY,my_BDY)  ) 

status = NF90_INQ_VARID(fidBDY,"Bathymetry",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_BDY")
status = NF90_INQ_VARID(fidBDY,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_BDY")
status = NF90_INQ_VARID(fidBDY,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_BDY")

status = NF90_GET_VAR(fidBDY,Bathymetry_ID,Bathymetry_BDY); call erreur(status,.TRUE.,"getvar_Bathymetry_BDY")
status = NF90_GET_VAR(fidBDY,nav_lat_ID,nav_lat_BDY);       call erreur(status,.TRUE.,"getvar_nav_lat_BDY")
status = NF90_GET_VAR(fidBDY,nav_lon_ID,nav_lon_BDY);       call erreur(status,.TRUE.,"getvar_nav_lon_BDY")
   
status = NF90_CLOSE(fidBDY); call erreur(status,.TRUE.,"close_coarse_bathy_file")

    
!===================================================================
! put COARSE grid in a npts-pts halo (+ transition in another halo)
npts=3  ! in nb of points on the regional grid

!=====================================================================
! 2- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
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
! 3- Extract variables on the REG grid :
!=================================================================================
    
ALLOCATE(  nav_lat_REG        (mx_REG,my_REG)  )
ALLOCATE(  nav_lon_REG        (mx_REG,my_REG)  )
if ( ln_isfcav) then
  ALLOCATE(  isf_draft_REG      (mx_REG,my_REG)  )
  ALLOCATE(  Bathymetry_isf_REG (mx_REG,my_REG)  )
endif
ALLOCATE(  Bathymetry_REG     (mx_REG,my_REG)  )
    
nav_lat_REG        (1:mx_REG,1:my_REG) = nav_lat_GLO        (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
nav_lon_REG        (1:mx_REG,1:my_REG) = nav_lon_GLO        (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
if ( ln_isfcav) then
  isf_draft_REG      (1:mx_REG,1:my_REG) = isf_draft_GLO      (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
  Bathymetry_isf_REG (1:mx_REG,1:my_REG) = Bathymetry_isf_GLO (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
endif
Bathymetry_REG     (1:mx_REG,1:my_REG) = Bathymetry_GLO     (imin_ORCA12:imax_ORCA12,jmin_ORCA12:jmax_ORCA12)
    
!---------------------------------------
! Adjust bathy along the edges of the REG grid :

write(*,*) 'Halo with bathy from the dataset used as lateral boundaries...'

!=== Interpolate the bathymetry of lateral boundary conditions over a npts-point halo :
do jREG=1,my_REG

    do iREG=1,npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

    do iREG=mx_REG-npts+1,mx_REG
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

enddo

do iREG=npts+1,mx_REG-npts

    do jREG=1,npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

    do jREG=my_REG-npts+1,my_REG
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

enddo

write(*,*) 'Smooth transition...'

!=== smooth transition from the Bathymetry of lateral boundaries to the interior bathymetry (over npts points again) :
do jREG=npts+1,my_REG-npts

    do iREG=npts+1,2*npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (2*npts+1-iREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
        &                               + (iREG-npts)     * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

    do iREG=mx_REG-2*npts+1,mx_REG-npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (mx_REG-npts+1-iREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
        &                               + (iREG-mx_REG+2*npts) * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

enddo

do iREG=2*npts+1,mx_REG-2*npts

    do jREG=npts+1,2*npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (2*npts+1-jREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
        &                               + (jREG-npts)     * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

    do jREG=my_REG-2*npts+1,my_REG-npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (my_REG-npts+1-jREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
        &                               + (jREG-my_REG+2*npts) * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
      endif
      if ( ln_isfcav) then
        isf_draft_REG      (iREG,jREG) = 0.e0
        Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
      endif
    enddo

enddo
   
!=================================================================================
! 4- Manual corrections for WED12 :
!=================================================================================

if ( TRIM(config) == 'WED12' ) then
   
    ! correction to avoid a closed cavity of 2x2x2 pts (identified after first mesh_mask creation)
    isf_draft_REG     (241:242,667:668) = 0.0 
    Bathymetry_isf_REG(241:242,667:668) = 0.0

    ! no isf along eastern boundary :
    isf_draft_REG     (1095:1122,668:703) = 0.0
    Bathymetry_isf_REG(1095:1122,668:703) = Bathymetry_REG(1095:1122,668:703)
    
    ! boxes to fill the Bellingshausen Sea : filled over [imin:imax,jmin:my]
    imin = (/   1 , 192 , 213 , 237 , 254 , 275 , 287 , 299 /)
    imax = (/ 191 , 212 , 236 , 253 , 274 , 286 , 298 , 325 /)
    jmin = (/ 494 , 807 , 835 , 862 , 876 , 894 , 899 , 903 /)
    write(*,*) 'Note for future bdy building:'
    write(*,*) '                    '
    write(*,*) '  ii_bdy_west(1)  = ', imax(8)+1
    write(*,*) '  j1_bdy_west(1)  = ', jmin(8)
    write(*,*) '  j2_bdy_west(1)  = ', my_REG-1 ! given that jmax(8)=my_REG-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(1) = ', imin(8)
    write(*,*) '  i2_bdy_north(1) = ', imax(8)
    write(*,*) '  jj_bdy_north(1) = ', jmin(8)-1
    write(*,*) '                    '
    write(*,*) '  i1_bdy_north(2) = ', imax(8)+2
    write(*,*) '  i2_bdy_north(2) = ', mx_REG-2 ! -2 because east bdy already contains mx_REG-1
    write(*,*) '  jj_bdy_north(2) = ', my_REG-1
    write(*,*) '                    '

    !----- put BDY along modified North-Western corner :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+1,imax(8)+npts
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
          endif
          if ( ln_isfcav) then
            isf_draft_REG      (iREG,jREG) = 0.e0
            Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          endif
      enddo
    enddo

    !- bdy_north(1) :
    do jREG=jmin(8)-npts, jmin(8)-1
      do iREG=imin(8),imax(8)
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
            &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
          endif
          if ( ln_isfcav) then
            isf_draft_REG      (iREG,jREG) = 0.e0
            Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          endif
      enddo
    enddo

    !----- smooth transition :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+npts+1,imax(8)+2*npts
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   (imax(8)+2*npts-iREG+1) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
            &                               + (iREG-imax(8)-npts)     * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
          endif
          if ( ln_isfcav) then
            isf_draft_REG      (iREG,jREG) = 0.e0
            Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          endif
      enddo
    enddo

    !- bdy_north(1) :
    do jREG=jmin(8)-2*npts,jmin(8)-npts-1
      do iREG=imin(8),imax(8)
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   (jmin(8)-npts-jREG)     * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
            &                               + (jREG-jmin(8)+2*npts+1) * Bathymetry_REG(iREG,jREG)                                                            ) / (npts+1)  
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
          endif
          if ( ln_isfcav) then
            isf_draft_REG      (iREG,jREG) = 0.e0
            Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          endif
      enddo
    enddo

    !- fill bellingshausen:
    do kk=1,8
      isf_draft_REG      (imin(kk):imax(kk),jmin(kk):my_REG) = 0.e0
      Bathymetry_isf_REG (imin(kk):imax(kk),jmin(kk):my_REG) = 0.e0
      Bathymetry_REG     (imin(kk):imax(kk),jmin(kk):my_REG) = 0.e0
    enddo

endif

!=================================================================================
! 5- Writing new REG bathymetry file :
!=================================================================================

write(*,*) 'Creating ', TRIM(file_bathy_out)

status = NF90_CREATE(TRIM(file_bathy_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')                     

status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")

status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID); call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID); call erreur(status,.TRUE.,"def_var_nav_lon_ID")
if ( ln_isfcav) then
  status = NF90_DEF_VAR(fidM,"isf_draft",NF90_FLOAT,(/dimID_x,dimID_y/),isf_draft_ID); call erreur(status,.TRUE.,"def_var_isf_draft_ID")
  status = NF90_DEF_VAR(fidM,"Bathymetry_isf",NF90_FLOAT,(/dimID_x,dimID_y/),Bathymetry_isf_ID); call erreur(status,.TRUE.,"def_var_Bathymetry_isf_ID")
endif
status = NF90_DEF_VAR(fidM,"Bathymetry",NF90_FLOAT,(/dimID_x,dimID_y/),Bathymetry_ID); call erreur(status,.TRUE.,"def_var_Bathymetry_ID")

status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"longname","Latitude")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"standard_name","latitude")
call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"longname","Longitude")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"standard_name","longitude")
call erreur(status,.TRUE.,"put_att_nav_lon_ID")
if ( ln_isfcav) then
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"long_name","ice-shelf draft")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"units","m")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,isf_draft_ID,"coordinates","nav_lat nav_lon")
  call erreur(status,.TRUE.,"put_att_isf_draft_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"long_name","bathymetry")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"units","m")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
  status = NF90_PUT_ATT(fidM,Bathymetry_isf_ID,"coordinates","nav_lat nav_lon")
  call erreur(status,.TRUE.,"put_att_Bathymetry_isf_ID")
endif
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"long_name","bathymetry with masked ice-shelf cavities")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"units","m")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")
status = NF90_PUT_ATT(fidM,Bathymetry_ID,"coordinates","nav_lat nav_lon")
call erreur(status,.TRUE.,"put_att_Bathymetry_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bathy.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL1")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"domain",TRIM(aaa))             ; call erreur(status,.TRUE.,"put_att_GLOBAL2")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,nav_lat_REG); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,nav_lon_REG); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( ln_isfcav) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_REG);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_REG); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_REG); call erreur(status,.TRUE.,"var_Bathymetry_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

DEALLOCATE( nav_lat_GLO, nav_lon_GLO, isf_draft_GLO, Bathymetry_isf_GLO, Bathymetry_GLO)
DEALLOCATE( nav_lat_BDY, nav_lon_BDY, Bathymetry_BDY)
DEALLOCATE( nav_lat_REG, nav_lon_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG)

write(*,*) '[done]'

end program modif


!=======================================================
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
    WRITE(*,*) 'ERROR:   ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'WHICH MEANS: ',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
