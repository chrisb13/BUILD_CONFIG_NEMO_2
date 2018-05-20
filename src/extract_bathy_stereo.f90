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
! 1- Read stereographic bathymetry and ice shelf draft
! 2- Read grid correspondance with GLO (i.e. extraction coordinates)
! 3- Read coarse bathymetry used for consistent bathymetry along boundaries
! 4- READ INTERPOLATION COEFFICIENTS FOR REGIONAL CONFIGURATION
! 5- Calculate bathy/isf draft on the REG grid
! 6- Writing new REG bathymetry file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_in_bathy_bdy, nn_isfcav,      &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, file_in_coord_bdy, ln_dateline, nn_perio
namelist /bathy_special/ file_spe_bathy, file_spe_isf_draft
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract,   &
&                                        nn_perio, nn_isfcav
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir, &
&                                        file_in_coord_bdy, file_spe_bathy, file_spe_isf_draft
LOGICAL                               :: ln_dateline

!-- STEREO variables :
INTEGER :: fidSTEREO1, fidSTEREO2, x_ID, y_ID, my_STEREO, mx_STEREO, bedrock_topography_ID, lat_ID, lon_ID, &
&          bathy_STEREO_ID, isf_draft_STEREO_ID, rq

REAL*8 :: a, e, lat_c, pm, lon_0, lon_STEREO, chi, m_c, t_c, t, x, y, res_STEREO, res_REG

REAL*4,ALLOCATABLE,DIMENSION(:) :: x_STEREO, y_STEREO

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: bathy_STEREO, isf_draft_STEREO, lat_STEREO, zlon_STEREO


!-- local variables :
INTEGER :: fidORCA12, fidBDY, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, isf_draft_ID, Bathymetry_isf_ID, Bathymetry_ID, &
&          my_ORCA12, mx_ORCA12,  my_BDY, mx_BDY,  my_REG, mx_REG, imin_ORCA12, imax_ORCA12, jmin_ORCA12, jmax_ORCA12, npts,          &
&          fidCOORDreg, fidCOORDpar, minlon, maxlon, minlat, maxlat, imin_STEREO, imax_STEREO, jmin_STEREO, jmax_STEREO, iREG, jREG,  &
&          iSTEREO, jSTEREO, iREGm1, iREGp1, jREGm1, jREGp1, kk, mx_tmp, my_tmp, i0, j0, rs, imin_STEREO2, imax_STEREO2,              &
&          jmin_STEREO2, jmax_STEREO2

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
REAL*8                                 :: eps, pi, rad2deg, deg2rad, aa
 
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
READ (UNIT=1, NML=bathy_special)
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

pi = ABS(ACOS(-1.d0))
deg2rad=pi/180.d0
rad2deg=180d0/pi

if ( nn_isfcav .eq. 0 ) then
  write(*,*) 'Seriously ?! Why using this dataset in a configuration with no ice shelves ??'
  write(*,*) 'Think about it, and come back when you change your mind...'
  write(*,*) '(or adapt extract_bathy_stereo.f90)'
  stop
endif

!=================================================================================
! 1- Read bathymetry and ice shelf draft on the stereographic grid
!=================================================================================

write(*,*) 'Reading bathy in ', TRIM(file_spe_bathy)

status = NF90_OPEN(TRIM(file_spe_bathy),0,fidSTEREO1); call erreur(status,.TRUE.,"Sart read STEREO") 

status = NF90_INQ_DIMID(fidSTEREO1,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidSTEREO1,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")

status = NF90_INQUIRE_DIMENSION(fidSTEREO1,dimID_y,len=my_STEREO); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidSTEREO1,dimID_x,len=mx_STEREO); call erreur(status,.TRUE.,"inq_dim_x")

ALLOCATE(  bathy_STEREO     (mx_STEREO,my_STEREO)  )
ALLOCATE(  isf_draft_STEREO (mx_STEREO,my_STEREO)  ) 
ALLOCATE(  x_STEREO (mx_STEREO)  ) 
ALLOCATE(  y_STEREO (my_STEREO)  ) 

status = NF90_INQ_VARID(fidSTEREO1,"bedrock",bathy_STEREO_ID);  call erreur(status,.TRUE.,"inq_bathy_STEREO_ID")
status = NF90_INQ_VARID(fidSTEREO1,"y",y_ID);                   call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidSTEREO1,"x",x_ID);                   call erreur(status,.TRUE.,"inq_lon_ID")

status = NF90_GET_VAR(fidSTEREO1,bathy_STEREO_ID,bathy_STEREO); call erreur(status,.TRUE.,"getvar_bathy_STEREO")
status = NF90_GET_VAR(fidSTEREO1,y_ID,y_STEREO);                call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidSTEREO1,x_ID,x_STEREO);                call erreur(status,.TRUE.,"getvar_lon")

status = NF90_CLOSE(fidSTEREO1); call erreur(status,.TRUE.,"End read bathy STEREO")

!-----

write(*,*) 'Converting x,y to lon,lat'

a     = 6378137.00000000
e     =       0.08181919
lat_c =     -71.00000000
lon_0 =       0.d0

write(*,*) '        Earth Radius                          = ', a
write(*,*) '        Earth misshapenness (excentricity)    = ', e
write(*,*) '        latitude of true scale in degrees     = ', lat_c
write(*,*) '        meridian in along the positive Y axis = ', lon_0 

!- if the standard parallel is in S.Hemi., switch signs.
if ( lat_c .lt. 0.0 ) then
  pm     = -1.d0    ! plus or minus, north lat. or south
  lat_c  = -lat_c
  lon_0  = -lon_0
else
  pm     = 1.d0
endif

!- convert to radians :
lat_c    = deg2rad * lat_c
lon_0    = deg2rad * lon_0

ALLOCATE( lat_STEREO(mx_STEREO,my_STEREO), zlon_STEREO(mx_STEREO,my_STEREO) )

do iSTEREO=1,mx_STEREO
write(*,*) iSTEREO
do jSTEREO=1,my_STEREO

  if ( pm .lt. 0.0 ) then
    x      = -x_STEREO(iSTEREO)
    y      = -y_STEREO(jSTEREO)
  else
    x      =  x_STEREO(iSTEREO)
    y      =  y_STEREO(jSTEREO)
  endif

  !- See Snyder for details.
  t_c  = tan(pi/4-lat_c/2) / ( (1-e*sin(lat_c)) / (1+e*sin(lat_c)))**(e/2)
  m_c  = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )
  t    = sqrt(x**2+y**2) * t_c / ( a * m_c )

  chi = 0.5*pi - 2 * atan(t) !- find lat with a series instead of iterating.

  lat_STEREO(iSTEREO,jSTEREO) = chi + ( (1./2.) * e**2 + (5./24.) * e**4 + ( 1./ 12.) * e**6 + (  13./   360.) * e**8 ) * sin(2*chi) &
  &                                 + (                  (7./48.) * e**4 + (29./240.) * e**6 + ( 811./ 11520.) * e**8 ) * sin(4*chi) &
  &                                 + (                                  + ( 7./120.) * e**6 + (  81./  1120.) * e**8 ) * sin(6*chi) &
  &                                 + (                                                      + (4279./161280.) * e**8 ) * sin(8*chi)

  lon_STEREO = lon_0 + atan2(x,-y)

  !- correct the signs and phasing :
  lat_STEREO(iSTEREO,jSTEREO) = pm * lat_STEREO(iSTEREO,jSTEREO)
  lon_STEREO                  = pm * lon_STEREO
  lon_STEREO                  = mod(lon_STEREO+pi,2*pi)-pi !- want longitude in the range -pi to pi

  !- convert back to degrees :
  lat_STEREO (iSTEREO,jSTEREO) = rad2deg * lat_STEREO(iSTEREO,jSTEREO)
  lon_STEREO                   = rad2deg * lon_STEREO

  zlon_STEREO(iSTEREO,jSTEREO) = lon_STEREO
  if ( ln_dateline ) then
    if( lon_STEREO .lt. 0.e0 ) then
      zlon_STEREO(iSTEREO,jSTEREO) = 360.e0 + lon_STEREO
    endif
  else
    if( lon_STEREO .gt. 180.e0 ) then
      zlon_STEREO(iSTEREO,jSTEREO) = lon_STEREO - 360.e0
    endif
  endif

enddo
enddo

!check :
write(*,*) 'min lon = ', MINVAL(zlon_STEREO)
write(*,*) 'max lon = ', MAXVAL(zlon_STEREO)
write(*,*) 'min lat = ', MINVAL( lat_STEREO)
write(*,*) 'max lat = ', MAXVAL( lat_STEREO)

!-----

write(*,*) 'Reading ice shelf draft in ', TRIM(file_spe_isf_draft)

status = NF90_OPEN(TRIM(file_spe_isf_draft),0,fidSTEREO2); call erreur(status,.TRUE.,"read isf_draft STEREO") 
status = NF90_INQ_VARID(fidSTEREO2,"ice_draft",isf_draft_STEREO_ID); call erreur(status,.TRUE.,"inq_isf_draft_STEREO_ID")
status = NF90_GET_VAR(fidSTEREO2,isf_draft_STEREO_ID,isf_draft_STEREO); call erreur(status,.TRUE.,"getvar_isf_draft_STEREO")
status = NF90_CLOSE(fidSTEREO2); call erreur(status,.TRUE.,"End read isf_draft STEREO")     

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
  
nn(:,:) = 0 
Bathymetry_isf_REG (:,:) = 0.e0
isf_draft_REG      (:,:) = 0.e0
Bathymetry_REG     (:,:) = 0.e0

res_REG = a * deg2rad * ABS( gphit_REG(1,2) - gphit_REG(1,1) ) 
res_STEREO = ABS( y_STEREO(2) - y_STEREO(1) )
rq = MIN( mx_STEREO , MIN( my_STEREO, 20 * CEILING( res_REG / res_STEREO ) ) )
write(*,*) ' '
write(*,*) 'Once a correspondance (iSTEREO,jSTEREO) on the stereographic grid is found for (iREG,jREG),'
write(*,*) ' the next stereographic point of the regional grid is searched in a square of'
write(*,*) ' radius ', rq, ' around (iSTEREO,jSTEREO).'
write(*,*) ' '

imin_STEREO2 = 1 
imax_STEREO2 = mx_STEREO
jmin_STEREO2 = 1
jmax_STEREO2 = my_STEREO

do iREG=1,mx_REG
!write(*,*) 'iREG = ', iREG
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

  ! restrict search for corresponding points on the stereographic grid (RAZ after a column)
  if ( ( jREG .eq. 1 ) .or. ( nn_perio .eq. 1 .and. ( iREG .lt. 11 .or. iREG .gt. mx_REG - 10 ) ) ) then
    imin_STEREO = 1
    imax_STEREO = mx_STEREO
    jmin_STEREO = 1
    jmax_STEREO = my_STEREO
  else
    imin_STEREO = MAX(         1, imin_STEREO2 )
    imax_STEREO = MIN( mx_STEREO, imax_STEREO2 )
    jmin_STEREO = MAX(         1, jmin_STEREO2 )
    jmax_STEREO = MIN( my_STEREO, jmax_STEREO2 )
  endif

  do iSTEREO=imin_STEREO,imax_STEREO
  do jSTEREO=jmin_STEREO,jmax_STEREO
 
    !-NB: we don't care too much about i=1 and i=mx (same for j) because it is masked anyway...    
    if (       zlon_STEREO(iSTEREO,jSTEREO) .ge. zglamt_REG(iREG,jREG) - 0.5*(zglamt_REG(iREG,jREG)-zglamt_REG(iREGm1,jREG)) &
    &    .and. zlon_STEREO(iSTEREO,jSTEREO) .lt. zglamt_REG(iREG,jREG) + 0.5*(zglamt_REG(iREGp1,jREG)-zglamt_REG(iREG,jREG)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .ge.  gphit_REG(iREG,jREG) - 0.5*( gphit_REG(iREG,jREG)- gphit_REG(iREG,jREGm1)) &
    &    .and.  lat_STEREO(iSTEREO,jSTEREO) .lt.  gphit_REG(iREG,jREG) + 0.5*( gphit_REG(iREG,jREGp1)- gphit_REG(iREG,jREG)) ) then

      if ( nn(iREG,jREG) .eq. 0 ) then
        imin_STEREO2 = mx_STEREO ! starting value for the MIN
        imax_STEREO2 = 1         ! starting value for the MAX
        jmin_STEREO2 = my_STEREO
        jmax_STEREO2 = 1
      endif

      Bathymetry_isf_REG (iREG,jREG) = Bathymetry_isf_REG (iREG,jREG) - bathy_STEREO     (iSTEREO,jSTEREO)
      isf_draft_REG      (iREG,jREG) = isf_draft_REG      (iREG,jREG) - isf_draft_STEREO (iSTEREO,jSTEREO)

      nn(iREG,jREG) = nn(iREG,jREG) + 1

      imin_STEREO2 = MIN( imin_STEREO2, iSTEREO - rq )
      imax_STEREO2 = MAX( imax_STEREO2, iSTEREO + rq )
      jmin_STEREO2 = MIN( jmin_STEREO2, jSTEREO - rq )
      jmax_STEREO2 = MAX( jmax_STEREO2, jSTEREO + rq )

    endif 

  enddo
  enddo

enddo
enddo

!-----

where ( nn(:,:) .ge. 1 )
  Bathymetry_isf_REG (:,:) = Bathymetry_isf_REG (:,:) / nn(:,:)
  isf_draft_REG      (:,:) = isf_draft_REG      (:,:) / nn(:,:)
elsewhere
  Bathymetry_isf_REG (:,:) = 99999999.9
  isf_draft_REG      (:,:) = 99999999.9
endwhere

!-----
! Interpolation at points with nn=0 :
! Only look for first neighbour
! needs to be adapted if many missing points

do iREG=1,mx_REG
do jREG=1,my_REG

  if ( nn(iREG,jREG) .eq. 0 ) then

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if    ( iREG+rs .gt. mx_REG ) then
          iREGp1 = 2+rs
        else
          iREGp1 = iREG + rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if    ( iREG+rs .gt. mx_REG ) then
          iREGp1 = mx_REG
          exit ! with nn(iREGp1,jREG)=0
        else
          iREGp1 = iREG + rs
        endif
      endif
      if ( nn(iREGp1,jREG) .ne. 0 ) exit ! with nn(iREGp1,jREG)>0
    enddo

    do rs=1,10
      if ( nn_perio .eq. 1 ) then !- periodic
        if ( iREG-rs .lt. 1     ) then
          iREGm1 = mx_REG-rs-1
        else
          iREGm1 = iREG - rs
        endif
      elseif ( nn_perio .eq. 0 ) then
        if ( iREG-rs .lt. 1     ) then
          iREGm1 = 1
          exit ! with nn(iREGm1,jREG)=0
        else
          iREGm1 = iREG - rs
        endif
      endif
      if ( nn(iREGm1,jREG) .gt. 0 ) exit ! with nn(iREGm1,jREG)>0
    enddo

    do rs=1,10
      if    ( jREG+rs .gt. my_REG ) then
        jREGp1 = my_REG
        exit ! with nn(iREG,jREGp1)=0
      else
        jREGp1 = jREG + rs
      endif
      if ( nn(iREG,jREGp1) .gt. 0 ) exit ! with nn(iREG,jREGp1)>0
    enddo

    do rs=1,10
      if    ( jREG-rs .lt. 1 ) then
        jREGm1 = 1
        exit ! with nn(iREG,jREGm1)=0
      else
        jREGm1 = jREG - rs
      endif
      if ( nn(iREG,jREGm1) .gt. 0 ) exit ! with nn(iREG,jREGm1)>0
    enddo

    Bathymetry_isf_REG(iREG,jREG) = (   nn(iREGm1,jREG  ) * (iREGp1-iREG) * Bathymetry_isf_REG(iREGm1,jREG  )   &
    &                                 + nn(iREGp1,jREG  ) * (iREG-iREGm1) * Bathymetry_isf_REG(iREGp1,jREG  )   &
    &                                 + nn(iREG  ,jREGm1) * (jREGp1-jREG) * Bathymetry_isf_REG(iREG  ,jREGm1)   &
    &                                 + nn(iREG  ,jREGp1) * (jREG-jREGm1) * Bathymetry_isf_REG(iREG  ,jREGp1) ) &
    &                             / (   nn(iREGm1,jREG  ) * (iREGp1-iREG)                                       &
    &                                 + nn(iREGp1,jREG  ) * (iREG-iREGm1)                                       &
    &                                 + nn(iREG  ,jREGm1) * (jREGp1-jREG)                                       &
    &                                 + nn(iREG  ,jREGp1) * (jREG-jREGm1)                                     )

    isf_draft_REG(iREG,jREG) = (   nn(iREGm1,jREG  ) * (iREGp1-iREG) * isf_draft_REG(iREGm1,jREG  )   &
    &                            + nn(iREGp1,jREG  ) * (iREG-iREGm1) * isf_draft_REG(iREGp1,jREG  )   &
    &                            + nn(iREG  ,jREGm1) * (jREGp1-jREG) * isf_draft_REG(iREG  ,jREGm1)   &
    &                            + nn(iREG  ,jREGp1) * (jREG-jREGm1) * isf_draft_REG(iREG  ,jREGp1) ) &
    &                        / (   nn(iREGm1,jREG  ) * (iREGp1-iREG)                                       &
    &                            + nn(iREGp1,jREG  ) * (iREG-iREGm1)                                       &
    &                            + nn(iREG  ,jREGm1) * (jREGp1-jREG)                                       &
    &                            + nn(iREG  ,jREGp1) * (jREG-jREGm1)                                     )

   endif

enddo
enddo

!---------

write(*,*) '>i=125 ', nn(125,20), Bathymetry_isf_REG(125,20), isf_draft_REG(125,20)
write(*,*) '>i=624 ', nn(624,20), Bathymetry_isf_REG(624,20), isf_draft_REG(624,20)

where( isf_draft_REG(:,:) .lt. 0.e0 )
  isf_draft_REG(:,:) = 0.e0
endwhere

where( Bathymetry_isf_REG(:,:) .lt. 0.e0 )
  Bathymetry_isf_REG(:,:) = 0.e0
endwhere

where( isf_draft_REG(:,:) .gt. 1.e0 .and. isf_draft_REG(:,:) .lt. 1.e4 )
  Bathymetry_REG(:,:) = 0.e0
elsewhere
  Bathymetry_REG(:,:) = Bathymetry_isf_REG(:,:)
endwhere

write(*,*) '#i=125 ', Bathymetry_REG(125,20), Bathymetry_isf_REG(125,20), isf_draft_REG(125,20)
write(*,*) '#i=624 ', Bathymetry_REG(624,20), Bathymetry_isf_REG(624,20), isf_draft_REG(624,20)

!---------------------------------------
! Adjust bathy along the edges of the REG grid :

write(*,*) 'Halo with bathy from the dataset used as lateral boundaries...'

!=== Interpolate the bathymetry of lateral boundary conditions over a npts-point halo :
do jREG=1,my_REG

    !-- Western BDY :
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

    !-- Eastern BDY :
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
! 6- Writing new REG bathymetry file :
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
if ( nn_isfcav .ge. 1 ) then
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
