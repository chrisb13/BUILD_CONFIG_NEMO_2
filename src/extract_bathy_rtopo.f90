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
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_in_bathy_bdy, nn_isfcav,      &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, file_in_coord_bdy, ln_dateline, nn_perio
namelist /rtopo/ file_rtopo_bathy, file_rtopo_isf_draft
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract,   &
&                                        nn_perio, nn_isfcav
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir, &
&                                        file_in_coord_bdy, file_rtopo_bathy, file_rtopo_isf_draft
LOGICAL                               :: ln_dateline

!-- RTOPO variables :
INTEGER :: fidRTOPO1, fidRTOPO2, dimID_latdim, dimID_londim, my_RTOPO, mx_RTOPO, bedrock_topography_ID, lat_ID, lon_ID, &
&          bathy_RTOPO_ID, isf_draft_RTOPO_ID

REAL*4,ALLOCATABLE,DIMENSION(:) :: lat_RTOPO, lon_RTOPO, zlon_RTOPO, tmp1, tmp2, tmp3, tmp4

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: bathy_RTOPO, isf_draft_RTOPO


!-- local variables :
INTEGER :: fidORCA12, fidBDY, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_ORCA12, mx_ORCA12,  my_BDY, mx_BDY,  my_REG, mx_REG, imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12, npts,          &
&          fidCOORDreg, fidCOORDpar, minlon, maxlon, minlat, maxlat, imin_RTOPO, imax_RTOPO, jmin_RTOPO, jmax_RTOPO, iREG, jREG,      &
&          iRTOPO, jRTOPO, iREGm1, iREGp1, jREGm1, jREGp1, kk, mx_tmp, my_tmp, i0, j0

CHARACTER(LEN=150) :: file_bathy_out, file_in_coord_REG

INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nn

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: gphit_REG, glamt_REG, zglamt_REG, isf_draft_REG, Bathymetry_isf_REG, Bathymetry_REG, &
&                                          Bathymetry_BDY, Bathymetry_isf_BDY, isf_draft_BDY

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
nn_isfcav         = 0

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
READ (UNIT=1, NML=rtopo)
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

if ( nn_isfcav .eq. 0 ) then
  write(*,*) 'Seriously ?! Why using RTopo in a configuration with no ice shelves ??'
  write(*,*) 'Think about it, and come back when you change your mind...'
  stop
endif

!=================================================================================
! 1- Read RTopo bathymetry and ice shelf draft
!=================================================================================

write(*,*) 'Reading RTopo bathy in ', TRIM(file_rtopo_bathy)

status = NF90_OPEN(TRIM(file_rtopo_bathy),0,fidRTOPO1); call erreur(status,.TRUE.,"Sart read RTopo") 

status = NF90_INQ_DIMID(fidRTOPO1,"latdim",dimID_latdim); call erreur(status,.TRUE.,"inq_dimID_latdim")
status = NF90_INQ_DIMID(fidRTOPO1,"londim",dimID_londim); call erreur(status,.TRUE.,"inq_dimID_londim")

status = NF90_INQUIRE_DIMENSION(fidRTOPO1,dimID_latdim,len=my_RTOPO); call erreur(status,.TRUE.,"inq_dim_latdim")
status = NF90_INQUIRE_DIMENSION(fidRTOPO1,dimID_londim,len=mx_RTOPO); call erreur(status,.TRUE.,"inq_dim_londim")

ALLOCATE(  bathy_RTOPO     (mx_RTOPO,my_RTOPO)  )
ALLOCATE(  isf_draft_RTOPO (mx_RTOPO,my_RTOPO)  ) 
ALLOCATE(  lat_RTOPO (my_RTOPO)  ) 
ALLOCATE(  lon_RTOPO (mx_RTOPO)  ) 
ALLOCATE( zlon_RTOPO (mx_RTOPO)  ) 

status = NF90_INQ_VARID(fidRTOPO1,"bedrock_topography",bathy_RTOPO_ID); call erreur(status,.TRUE.,"inq_bathy_RTOPO_ID")
status = NF90_INQ_VARID(fidRTOPO1,"lat",lat_ID);                        call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidRTOPO1,"lon",lon_ID);                        call erreur(status,.TRUE.,"inq_lon_ID")

status = NF90_GET_VAR(fidRTOPO1,bathy_RTOPO_ID,bathy_RTOPO); call erreur(status,.TRUE.,"getvar_bathy_RTOPO")
status = NF90_GET_VAR(fidRTOPO1,lat_ID,lat_RTOPO);           call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidRTOPO1,lon_ID,lon_RTOPO);           call erreur(status,.TRUE.,"getvar_lon")

status = NF90_CLOSE(fidRTOPO1); call erreur(status,.TRUE.,"End read bathy RTopo")     

zlon_RTOPO(:) = lon_RTOPO(:)
if ( ln_dateline ) then
  where( lon_RTOPO(:) .lt. 0.e0 )
    zlon_RTOPO(:) = 360.e0 + lon_RTOPO(:)
  endwhere
else
  where( lon_RTOPO(:) .gt. 180.e0 )
    zlon_RTOPO(:) = lon_RTOPO(:) - 360.e0
  endwhere
endif

!-----

write(*,*) 'Reading RTopo ice shelf draft in ', TRIM(file_rtopo_isf_draft)

status = NF90_OPEN(TRIM(file_rtopo_isf_draft),0,fidRTOPO2); call erreur(status,.TRUE.,"read isf_draft RTopo") 
status = NF90_INQ_VARID(fidRTOPO2,"ice_base_topography",isf_draft_RTOPO_ID); call erreur(status,.TRUE.,"inq_isf_draft_RTOPO_ID")
status = NF90_GET_VAR(fidRTOPO2,isf_draft_RTOPO_ID,isf_draft_RTOPO); call erreur(status,.TRUE.,"getvar_isf_draft_RTOPO")
status = NF90_CLOSE(fidRTOPO2); call erreur(status,.TRUE.,"End read isf_draft RTopo")     

!=================================================================================
! 2- Read grid correspondance with GLO (i.e. extraction coordinates)
!=================================================================================

write(*,*) 'Reading REGIONAL lon,lat in ', TRIM(file_in_coord_REG)

status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidCOORDreg); call erreur(status,.TRUE.,"read coord input")

status = NF90_INQ_DIMID(fidCOORDreg,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidCOORDreg,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
status = NF90_INQUIRE_DIMENSION(fidCOORDreg,dimID_y,len=my_REG); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidCOORDreg,dimID_x,len=mx_REG); call erreur(status,.TRUE.,"inq_dim_x_REG")

ALLOCATE(  gphit_REG (mx_REG,my_REG)  )
ALLOCATE(  glamt_REG (mx_REG,my_REG)  )
ALLOCATE( zglamt_REG (mx_REG,my_REG)  )

status = NF90_INQ_VARID(fidCOORDreg,"gphit",nav_lat_ID); call erreur(status,.TRUE.,"inq_gphit_REG_ID")
status = NF90_INQ_VARID(fidCOORDreg,"glamt",nav_lon_ID); call erreur(status,.TRUE.,"inq_glamt_REG_ID")

status = NF90_GET_VAR(fidCOORDreg,nav_lat_ID,gphit_REG); call erreur(status,.TRUE.,"getvar_gphit_REG")
status = NF90_GET_VAR(fidCOORDreg,nav_lon_ID,glamt_REG); call erreur(status,.TRUE.,"getvar_glamt_REG")

!status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "imin_extraction", imin_ORCA12); call erreur(status,.TRUE.,"read att1")
!status = NF90_GET_ATT(fidCOORDreg, NF90_GLOBAL, "jmin_extraction", jmin_ORCA12); call erreur(status,.TRUE.,"read att2")

status = NF90_CLOSE(fidCOORDreg)                       ; call erreur(status,.TRUE.,"end read fidCOORDreg")

zglamt_REG(:,:)=glamt_REG(:,:)
if ( ln_dateline ) then
  where( glamt_REG(:,:) .lt. 0.e0 )
    zglamt_REG(:,:)=glamt_REG(:,:)+360.e0
  endwhere
else
  where( glamt_REG(:,:) .gt. 180.e0 )
    zglamt_REG(:,:)=glamt_REG(:,:)-360.e0
  endwhere
endif

!- lat/lon REG boundaries to limit loops on RTOPO's lon/lat
minlon = MINVAL(zglamt_REG-0.5)
maxlon = MAXVAL(zglamt_REG+0.5)
minlat = MINVAL( gphit_REG-0.5)
maxlat = MAXVAL( gphit_REG+0.5)

ALLOCATE( tmp1(mx_RTOPO), tmp2(mx_RTOPO) )
ALLOCATE( tmp3(my_RTOPO), tmp4(my_RTOPO) )
do iRTOPO=1,mx_RTOPO
  tmp1(iRTOPO) = abs( zlon_RTOPO(iRTOPO) - minlon )
  tmp2(iRTOPO) = abs( zlon_RTOPO(iRTOPO) - maxlon )
enddo
do jRTOPO=1,my_RTOPO
  tmp3(jRTOPO) = abs(  lat_RTOPO(jRTOPO) - minlat )
  tmp4(jRTOPO) = abs(  lat_RTOPO(jRTOPO) - maxlat )
enddo
imin_RTOPO = MINLOC( tmp1(:), 1 ) - 1
imax_RTOPO = MINLOC( tmp2(:), 1 ) + 1
jmin_RTOPO = MINLOC( tmp3(:), 1 ) - 1
jmax_RTOPO = MINLOC( tmp4(:), 1 ) + 1
DEALLOCATE( lon_RTOPO, tmp1, tmp2, tmp3, tmp4 )

!=================================================================================
! 3- Read coarse bathymetry used for consistent bathymetry along boundaries
!=================================================================================

write(*,*) 'Reading coarse bathymetry for consistent boundaries: ', TRIM(file_in_bathy_bdy)

status = NF90_OPEN(TRIM(file_in_bathy_bdy),0,fidBDY); call erreur(status,.TRUE.,"read_coarse_bathymetry") 

status = NF90_INQ_DIMID(fidBDY,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y_BDY")
status = NF90_INQ_DIMID(fidBDY,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x_BDY")

status = NF90_INQUIRE_DIMENSION(fidBDY,dimID_y,len=my_BDY); call erreur(status,.TRUE.,"inq_dim_y_BDY")
status = NF90_INQUIRE_DIMENSION(fidBDY,dimID_x,len=mx_BDY); call erreur(status,.TRUE.,"inq_dim_x_BDY")

ALLOCATE(  Bathymetry_BDY (mx_BDY,my_BDY)  ) 
if ( nn_isfcav .eq. 2 ) then
  ALLOCATE(  Bathymetry_isf_BDY (mx_BDY,my_BDY)  )
  ALLOCATE(  isf_draft_BDY (mx_BDY,my_BDY)  )
endif

status = NF90_INQ_VARID(fidBDY,"nav_lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"lat",nav_lat_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"latitude",nav_lat_ID)
call erreur(status,.TRUE.,"inq_nav_lat_ID_BDY")
status = NF90_INQ_VARID(fidBDY,"nav_lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"lon",nav_lon_ID)
if ( status .ne. 0 ) status = NF90_INQ_VARID(fidBDY,"longitude",nav_lon_ID)
call erreur(status,.TRUE.,"inq_nav_lon_ID_BDY")
status = NF90_INQ_VARID(fidBDY,"Bathymetry",Bathymetry_ID)
call erreur(status,.TRUE.,"inq_Bathymetry_ID_BDY")
if ( nn_isfcav .eq. 2 ) then
  status = NF90_INQ_VARID(fidBDY,"Bathymetry_isf",Bathymetry_isf_ID)
  call erreur(status,.TRUE.,"inq_Bathymetry_isf_ID_BDY")
  status = NF90_INQ_VARID(fidBDY,"isf_draft",isf_draft_ID)
  call erreur(status,.TRUE.,"inq_isf_draft_ID_BDY")
endif

status = NF90_GET_VAR(fidBDY,Bathymetry_ID,Bathymetry_BDY);         call erreur(status,.TRUE.,"getvar_Bathymetry_BDY")
if ( nn_isfcav .eq. 2 ) then
  status = NF90_GET_VAR(fidBDY,Bathymetry_isf_ID,Bathymetry_isf_BDY); call erreur(status,.TRUE.,"getvar_Bathymetry_isf_BDY")
  status = NF90_GET_VAR(fidBDY,isf_draft_ID,isf_draft_BDY);           call erreur(status,.TRUE.,"getvar_isf_draft_BDY")
endif   

status = NF90_CLOSE(fidBDY); call erreur(status,.TRUE.,"close_coarse_bathy_file")

!===================================================================
! put COARSE grid in a npts-pts halo (+ transition in another halo)
npts=3  ! in nb of points on the regional grid

!=====================================================================
! 4- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
!    NB: must be consistent with the grid of file_in_bathy_bdy
!=====================================================================

write(*,*) 'Reading interpolation weights in ', TRIM(file_coeff)

status = NF90_OPEN(TRIM(file_coeff),0,fidcoeff) ; call erreur(status,.TRUE.,"read weights") 
                                                  
status = NF90_INQ_DIMID(fidcoeff,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidcoeff,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
                                                     
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_y,len=my_tmp); call erreur(status,.TRUE.,"inq_dim_y_REG")
status = NF90_INQUIRE_DIMENSION(fidcoeff,dimID_x,len=mx_tmp); call erreur(status,.TRUE.,"inq_dim_x_REG")
if ( mx_REG .ne. mx_tmp .or. my_REG .ne. my_tmp ) then
  write(*,*) '~!@#$%^* ERROR: mismatch between weights and regional coordinate files >>>>>>>>> stop !!'
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
! 5- Calculate bathy/isf draft on the REG grid :
!=================================================================================

write(*,*) 'Calculating bathymetry on the regional grid...'
    
ALLOCATE(  nn                 (mx_REG,my_REG)  )
ALLOCATE(  isf_draft_REG      (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_isf_REG (mx_REG,my_REG)  )
ALLOCATE(  Bathymetry_REG     (mx_REG,my_REG)  )
   
Bathymetry_isf_REG (:,:) = 0.e0
isf_draft_REG      (:,:) = 0.e0
Bathymetry_REG     (:,:) = 0.e0

do iREG=1,mx_REG
write(*,*) 'iREG = ', iREG
do jREG=1,my_REG

  if ( nn_perio .eq. 1 ) then !- periodic
    if    ( iREG+1 .gt. mx_REG ) then
      iREGp1 = 3
      iREGm1 = iREG - 1
    elseif ( iREG-1 .lt. 1     ) then
      iREGp1 = iREG + 1
      iREGm1 = mx_REG-2
    else
      iREGp1 = iREG + 1
      iREGm1 = iREG - 1
    endif
  elseif ( nn_perio .eq. 0 ) then
    if    ( iREG+1 .gt. mx_REG ) then
      iREGp1 = iREG
      iREGm1 = iREG - 1
    elseif ( iREG-1 .lt. 1     ) then
      iREGp1 = iREG + 1
      iREGm1 = iREG
    else
      iREGp1 = iREG + 1
      iREGm1 = iREG - 1
    endif
  else
    write(*,*) '~!@#$%^* nn_perio must be either 0 or 1 >>>>> stop !!'
    stop
  endif
  !-
  if    ( jREG+1 .gt. my_REG ) then
    jREGp1 = jREG
    jREGm1 = jREG - 1
  elseif ( iREG-1 .lt. 1     ) then
    jREGp1 = jREG + 1
    jREGm1 = jREG
  else
    jREGp1 = jREG + 1
    jREGm1 = jREG - 1
  endif
  !-
  iREGm1 = MAX( MIN( iREGm1, mx_REG ), 1 )
  iREGp1 = MAX( MIN( iREGp1, mx_REG ), 1 )
  jREGm1 = MAX( MIN( jREGm1, my_REG ), 1 )
  jREGp1 = MAX( MIN( jREGp1, my_REG ), 1 )
  !-
  
  do iRTOPO=imin_RTOPO,imax_RTOPO
  do jRTOPO=jmin_RTOPO,jmax_RTOPO
 
    !-NB: we don't care too much about i=1 and i=mx (same for j) because it is masked anyway...    
    if (       zlon_RTOPO(iRTOPO) .ge. zglamt_REG(iREG,jREG) - 0.5*(zglamt_REG(iREG,jREG)-zglamt_REG(iREGm1,jREG)) &
    &    .and. zlon_RTOPO(iRTOPO) .lt. zglamt_REG(iREG,jREG) + 0.5*(zglamt_REG(iREGp1,jREG)-zglamt_REG(iREG,jREG)) &
    &    .and.  lat_RTOPO(jRTOPO) .ge.  gphit_REG(iREG,jREG) - 0.5*( gphit_REG(iREG,jREG)- gphit_REG(iREG,jREGm1)) &
    &    .and.  lat_RTOPO(jRTOPO) .lt.  gphit_REG(iREG,jREG) + 0.5*( gphit_REG(iREG,jREGp1)- gphit_REG(iREG,jREG)) ) then

      Bathymetry_isf_REG (iREG,jREG) = Bathymetry_isf_REG (iREG,jREG) - MIN(0.0, bathy_RTOPO     (iRTOPO,jRTOPO))
      isf_draft_REG      (iREG,jREG) = isf_draft_REG      (iREG,jREG) - MIN(0.0, isf_draft_RTOPO (iRTOPO,jRTOPO))

      nn(iREG,jREG) = nn(iREG,jREG) + 1

    endif 

  enddo
  enddo
  !-
  where ( nn(:,:) .ge. 1 )
    Bathymetry_isf_REG (:,:) = Bathymetry_isf_REG (:,:) / nn(:,:)
    isf_draft_REG      (:,:) = isf_draft_REG      (:,:) / nn(:,:)
  elsewhere
    Bathymetry_isf_REG (:,:) = 99999999.
    isf_draft_REG      (:,:) = 99999999.
  endwhere
  !- 
  where( isf_draft_REG(:,:) .gt. 1.e0 )
    Bathymetry_REG(:,:) = 0.e0
  elsewhere
    Bathymetry_REG(:,:) = Bathymetry_isf_REG(:,:)
  endwhere

enddo
enddo
 
!---------------------------------------
! Adjust bathy along the edges of the REG grid :

write(*,*) 'Halo with bathy from the dataset used as lateral boundaries...'

!=== Interpolate the bathymetry of lateral boundary conditions over a npts-point halo :
do jREG=1,my_REG

    !-- Eastern BDY :
    do iREG=1,npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        if ( nn_isfcav .eq. 1 ) then
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          isf_draft_REG      (iREG,jREG) = 0.e0 
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

    !-- Western BDY :
    aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
    do iREG=mx_REG-npts+1,mx_REG
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        if ( nn_isfcav .eq. 1 ) then
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          isf_draft_REG      (iREG,jREG) = 0.e0 
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

enddo

do iREG=npts+1,mx_REG-npts

    !- Southern BDY :
    do jREG=1,npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        if ( nn_isfcav .eq. 1 ) then
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          isf_draft_REG      (iREG,jREG) = 0.e0 
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

    !- Northern BDY :
    do jREG=my_REG-npts+1,my_REG
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
        &                               + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        if ( nn_isfcav .eq. 1 ) then
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          isf_draft_REG      (iREG,jREG) = 0.e0 
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
          isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
          &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

enddo

write(*,*) 'Smooth transition...'

!=== smooth transition from the Bathymetry of lateral boundaries to the interior bathymetry (over npts points again) :
do jREG=npts+1,my_REG-npts

    !-- Western BDY :
    do iREG=npts+1,2*npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (2*npts+1-iREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)                   &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa            &
        &                               + (iREG-npts)     * Bathymetry_REG(iREG,jREG)                                                         ) / (npts+1)
        if ( nn_isfcav .eq. 1 ) then
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
          isf_draft_REG      (iREG,jREG) = 0.e0
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   (2*npts+1-iREG) * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa    &
          &                                   + (iREG-npts)     * Bathymetry_isf_REG(iREG,jREG)                                                         ) / (npts+1)  
          isf_draft_REG (iREG,jREG) =  (   (2*npts+1-iREG) * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa    &
          &                              + (iREG-npts)     * isf_draft_REG(iREG,jREG)                                                         ) / (npts+1)  
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

    !-- Eastern BDY :
    do iREG=mx_REG-2*npts+1,mx_REG-npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (iREG-mx_REG+2*npts) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa   &
        &                               + (mx_REG-npts+1-iREG) * Bathymetry_REG(iREG,jREG)                                                        ) / (npts+1) 
        if ( nn_isfcav .eq. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   (iREG-mx_REG+2*npts) * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
          &                                                              + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
          &                                                              + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
          &                                                              + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
          &                                   + (mx_REG-npts+1-iREG) * Bathymetry_isf_REG(iREG,jREG)                                                        ) / (npts+1)
          isf_draft_REG (iREG,jREG) =  (   (iREG-mx_REG+2*npts) * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
          &                                                         + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
          &                                                         + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
          &                                                         + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
          &                              + (mx_REG-npts+1-iREG) * isf_draft_REG(iREG,jREG)                                                        ) / (npts+1)
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

enddo

do iREG=2*npts+1,mx_REG-2*npts

    !-- Southern BDY :
    do jREG=npts+1,2*npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (2*npts+1-jREG) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)           &
        &                                                     + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)           &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)           &
        &                                                     + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa    &
        &                               + (jREG-npts)     * Bathymetry_REG(iREG,jREG)                                                        ) / (npts+1)  
        if ( nn_isfcav .eq. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   (2*npts+1-jREG) * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)           &
          &                                                         + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa    &
          &                                   + (jREG-npts)     * Bathymetry_isf_REG(iREG,jREG)                                                        ) / (npts+1)  
          isf_draft_REG (iREG,jREG) =  (   (2*npts+1-jREG) * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)           &
          &                                                    + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa    &
          &                              + (jREG-npts)     * isf_draft_REG(iREG,jREG)                                                        ) / (npts+1)
        endif
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

    !-- Northern BDY :
    do jREG=my_REG-2*npts+1,my_REG-npts
      aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
      if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
        Bathymetry_REG (iREG,jREG) =  (   (jREG-my_REG+2*npts) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)          &
        &                                                          + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa   &
        &                               + (my_REG-npts+1-jREG) * Bathymetry_REG(iREG,jREG)                                                        ) / (npts+1)  
        if ( nn_isfcav .eq. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        elseif ( nn_isfcav .eq. 2 ) then
          Bathymetry_isf_REG (iREG,jREG) =  (   (jREG-my_REG+2*npts) * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)         &
          &                                                              + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)         &
          &                                                              + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)         &
          &                                                              + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa  &
          &                                   + (my_REG-npts+1-jREG) * Bathymetry_isf_REG(iREG,jREG)                                                        ) / (npts+1)  
          isf_draft_REG (iREG,jREG) =  (   (jREG-my_REG+2*npts) * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)          &
          &                                                         + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)          &
          &                                                         + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)          &
          &                                                         + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa   &
          &                              + (my_REG-npts+1-jREG) * isf_draft_REG(iREG,jREG)                                                        ) / (npts+1) 
        endif 
      else
        Bathymetry_REG (iREG,jREG) = 0.e0
        if ( nn_isfcav .ge. 1 ) then
          isf_draft_REG      (iREG,jREG) = 0.e0
          Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
        endif
      endif
    enddo

enddo
   
!=================================================================================
! 4- Manual corrections for WED12 :
!=================================================================================

if ( TRIM(config) == 'WED12' ) then
   
    ! Note that you DO NOT have to change the following in you change the domain size
    ! through modifications of nn_imin_extract, nn_imax_extract, ... in the namelist 

    ! To keep the boxes at the same position:
    i0 = imin_ORCA12 - 2464
    j0 = jmin_ORCA12 -  151

    ! correction to avoid a closed cavity of 2x2x2 pts (identified after first mesh_mask creation)
    isf_draft_REG     (i0+241:i0+242,j0+667:j0+668) = 0.0 
    Bathymetry_isf_REG(i0+241:i0+242,j0+667:j0+668) = 0.0

    ! no isf along eastern boundary (adapt manually to adjust more accurately) :
    isf_draft_REG     (i0+1095:mx_REG,j0+668:j0+703) = 0.0
    Bathymetry_isf_REG(i0+1095:mx_REG,j0+668:j0+703) = Bathymetry_REG(i0+1095:mx_REG,j0+668:j0+703)
    
    ! boxes to fill the Bellingshausen Sea : filled over [imin:imax,jmin:my]
    imin = (/   1+i0 , 192+i0 , 213+i0 , 237+i0 , 254+i0 , 275+i0 , 287+i0 , 299+i0 /)
    imax = (/ 191+i0 , 212+i0 , 236+i0 , 253+i0 , 274+i0 , 286+i0 , 298+i0 , 325+i0 /)
    jmin = (/ 494+j0 , 807+j0 , 835+j0 , 862+j0 , 876+j0 , 894+j0 , 899+j0 , 903+j0 /)
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
            if ( nn_isfcav .eq. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            elseif ( nn_isfcav .eq. 2 ) then
              Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
              isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
            endif
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
            if ( nn_isfcav .ge. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            endif
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
            if ( nn_isfcav .eq. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            elseif ( nn_isfcav .eq. 2 ) then
              Bathymetry_isf_REG (iREG,jREG) =  (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
              &                                   + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
              isf_draft_REG (iREG,jREG) =  (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)   &
              &                              + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa
            endif 
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
            if ( nn_isfcav .ge. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            endif
          endif
      enddo
    enddo

    !----- smooth transition :
    !- bdy_west(1) :
    do jREG=jmin(8),my_REG
      do iREG=imax(8)+npts+1,imax(8)+2*npts
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   (iREG-imax(8)-npts)     * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
            &                               + (imax(8)+2*npts-iREG+1) * Bathymetry_REG(iREG,jREG)                                                        ) / (npts+1)
            if ( nn_isfcav .eq. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            elseif ( nn_isfcav .eq. 2 ) then 
              Bathymetry_isf_REG (iREG,jREG) =  (   (iREG-imax(8)-npts)     * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
              &                                   + (imax(8)+2*npts-iREG+1) * Bathymetry_isf_REG(iREG,jREG)                                                        ) / (npts+1)
              isf_draft_REG (iREG,jREG) =  (   (iREG-imax(8)-npts)     * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
              &                              + (imax(8)+2*npts-iREG+1) * isf_draft_REG(iREG,jREG)                                                        ) / (npts+1)
            endif
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
            if ( nn_isfcav .ge. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            endif
          endif
      enddo
    enddo

    !- bdy_north(1) :
    do jREG=jmin(8)-2*npts,jmin(8)-npts-1
      do iREG=imin(8),imax(8)
          aa = wgNWt(iREG,jREG) + wgSWt(iREG,jREG) + wgNEt(iREG,jREG) + wgSEt(iREG,jREG)
          if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
            Bathymetry_REG (iREG,jREG) =  (   (jREG-jmin(8)+2*npts+1) * (   Bathymetry_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
            &                                                             + Bathymetry_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
            &                               + (jmin(8)-npts-jREG)     * Bathymetry_REG(iREG,jREG)                                                        ) / (npts+1)  
            if ( nn_isfcav .eq. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            elseif ( nn_isfcav .eq. 2 ) then
              Bathymetry_isf_REG (iREG,jREG) =  (   (jREG-jmin(8)+2*npts+1) * (   Bathymetry_isf_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
              &                                                                 + Bathymetry_isf_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
              &                                   + (jmin(8)-npts-jREG)     * Bathymetry_isf_REG(iREG,jREG)                                                        ) / (npts+1)  
              isf_draft_REG (iREG,jREG) =  (   (jREG-jmin(8)+2*npts+1) * (   isf_draft_BDY( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)        &
              &                                                            + isf_draft_BDY( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG) ) / aa &
              &                              + (jmin(8)-npts-jREG)     * isf_draft_REG(iREG,jREG)                                                        ) / (npts+1)  
            endif
          else
            Bathymetry_REG (iREG,jREG) = 0.e0
            if ( nn_isfcav .ge. 1 ) then
              isf_draft_REG      (iREG,jREG) = 0.e0
              Bathymetry_isf_REG (iREG,jREG) = Bathymetry_REG (iREG,jREG)
            endif
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
if ( nn_isfcav .ge. 1 ) then
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
if ( nn_isfcav .ge. 1 ) then
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

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_bathy_rtopo.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
call erreur(status,.TRUE.,"put_att_GLOBAL1")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin_extraction",imin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL3")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax_extraction",imax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL4")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin_extraction",jmin_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL5")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax_extraction",jmax_ORCA12)  ; call erreur(status,.TRUE.,"put_att_GLOBAL6")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,nav_lat_ID,gphit_REG); call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,glamt_REG); call erreur(status,.TRUE.,"var_nav_lon_ID")
if ( nn_isfcav .eq. 2 ) then
  status = NF90_PUT_VAR(fidM,isf_draft_ID,isf_draft_REG);           call erreur(status,.TRUE.,"var_isf_draft_ID")
  status = NF90_PUT_VAR(fidM,Bathymetry_isf_ID,Bathymetry_isf_REG); call erreur(status,.TRUE.,"var_Bathymetry_isf_ID")
endif
status = NF90_PUT_VAR(fidM,Bathymetry_ID,Bathymetry_REG); call erreur(status,.TRUE.,"var_Bathymetry_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

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
