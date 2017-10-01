program buildweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! by N. Jourdain, on 23-MAY-2013, at CCRC-UNSW, Sydney                    !
!                                                                         !
! build coefficients to inpterpole a global simulation on a regional grid !
! (needed to build the initial state of regional config)                  !
!                                                                         !
! 0- Initializations                                                      !
! 2- Read the grid coordinates for the REGIONAL domain                    !
! 3- CALCULATE WEIGHTS TO INTERPOLATE GLOBAL FILEDS ONTO THE REGIONAL GRID!
! 4- WRITE EVERYTHING NEEDED FOR INTERPOLATION IN A NETCDF FILE           !
!                                                                         !
! hisotry: - Jul. 2017: adaptation to BUILD_CONFIG_NEMO_2                 !
!                                                                         !
! WARNING : for the moment, only work :                                   !
!           - if the vertical grid in global input domain is the same as  !
!             in the regional output domain                               !
!           - if the regional grid is finer or no more than 3-4 times     !
!             coarser than the global ocean grid                          !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf

IMPLICIT NONE

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, file_in_bathy_bdy, nn_isfcav,     &
& nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, file_in_coord_bdy, ln_dateline, nn_perio
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract,  &
&                                        nn_isfcav, nn_perio
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, config_dir, file_in_coord_bdy
LOGICAL                               :: ln_dateline

!--
INTEGER                                  :: fidhgr, dimID_y, dimID_x, my_GLO, mx_GLO, gphif_ID, gphiv_ID, &
&                                           gphiu_ID, gphit_ID, glamf_ID, glamv_ID, glamu_ID, glamt_ID, &
&                                           iGLO, jGLO
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: gphitGLO, glamtGLO, zglamtGLO,  &
&                                           gphifGLO, glamfGLO, zglamfGLO,  & 
&                                           gphiuGLO, glamuGLO, zglamuGLO,  & 
&                                           gphivGLO, glamvGLO, zglamvGLO
REAL*4                                   :: anginX, anginY, ra, rad

!--
INTEGER                                  :: fidgridout, mtout, jpkout, my_REG, mx_REG, jkout, jREG, iREG,    &
&                                           fmask_ID, dimID_z, tmask_ID, umask_ID, vmask_ID,                 &
&                                           ffout_ID, status, glamvREG_ID, glamtREG_ID, glamuREG_ID,         & 
&                                           gphifREG_ID, gphivREG_ID, gphiuREG_ID, gphitREG_ID, glamfREG_ID
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: gphitREG, glamtREG, zglamtREG
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: gphifREG, glamfREG, zglamfREG
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: gphiuREG, glamuREG, zglamuREG
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: gphivREG, glamvREG, zglamvREG
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: angoutXt, angoutYt, angleXt, angleYt
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: angoutXf, angoutYf, angleXf, angleYf
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: angoutXu, angoutYu, angleXu, angleYu
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: angoutXv, angoutYv, angleXv, angleYv

!-- INTERPOLATION COEFFICIENTS
INTEGER                                  :: ziWt_ID, ziEt_ID, zjSt_ID, zjNt_ID, zkUt_ID, wgNEt_ID, wgSEt_ID, &
&                                           wgSWt_ID, wgNWt_ID, angleXt_ID, angleYt_ID, fidcoeff
INTEGER                                  :: ziWf_ID, ziEf_ID, zjSf_ID, zjNf_ID, zkUf_ID, wgNEf_ID, wgSEf_ID, &
&                                           wgSWf_ID, wgNWf_ID, angleXf_ID, angleYf_ID
INTEGER                                  :: ziWu_ID, ziEu_ID, zjSu_ID, zjNu_ID, zkUu_ID, wgNEu_ID, wgSEu_ID, &
&                                           wgSWu_ID, wgNWu_ID, angleXu_ID, angleYu_ID
INTEGER                                  :: ziWv_ID, ziEv_ID, zjSv_ID, zjNv_ID, zkUv_ID, wgNEv_ID, wgSEv_ID, &
&                                           wgSWv_ID, wgNWv_ID, angleXv_ID, angleYv_ID
REAL*4                                   :: dist, distmin, distSW, distNW, distSE, distNE
INTEGER,ALLOCATABLE,DIMENSION(:,:)       :: ziWt, ziEt, zjSt, zjNt, zzit, zzjt, zzitm1, zzjtm1, zzitp1, zzjtp1
INTEGER,ALLOCATABLE,DIMENSION(:,:)       :: ziWf, ziEf, zjSf, zjNf, zzif, zzjf, zzifm1, zzjfm1, zzifp1, zzjfp1
INTEGER,ALLOCATABLE,DIMENSION(:,:)       :: ziWu, ziEu, zjSu, zjNu, zziu, zzju, zzium1, zzjum1, zziup1, zzjup1
INTEGER,ALLOCATABLE,DIMENSION(:,:)       :: ziWv, ziEv, zjSv, zjNv, zziv, zzjv, zzivm1, zzjvm1, zzivp1, zzjvp1
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: wgNEt, wgSEt, wgSWt, wgNWt
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: wgNEf, wgSEf, wgSWf, wgNWf
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: wgNEu, wgSEu, wgSWu, wgNWu
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: wgNEv, wgSEv, wgSWv, wgNWv
CHARACTER(LEN=150)                       :: file_coeff, file_in_coord_REG

!=================================================================================
! 0- Initializations 
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'
file_in_bathy_bdy = 'not_used'
nn_isfcav         = 0
ln_dateline       = .false.
nn_perio          = 0

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
CLOSE(1)

! name of regional coordinates file (output file) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

write(file_coeff,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coeff_linear_',a,'.nc')

!- earth radius (meter), must agree with NEMO file phycst.F90
ra = 6371229.0

!- deg to rad conversion (same as phycst.F90)
rad = 3.14159265358979323846264338327 / 180.000000000000000000000000000

!==================================================================================
! 1- Read the grid coordinate for the GLOBAL domain used as BDY conditions
!==================================================================================

write(*,*) 'Reading ', TRIM(file_in_coord_bdy)

status = NF90_OPEN(TRIM(file_in_coord_bdy),0,fidhgr)          
call erreur(status,.TRUE.,"read_horizontal_grid") 
                                       
status = NF90_INQ_DIMID(fidhgr,"y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidhgr,"x",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
                                           
status = NF90_INQUIRE_DIMENSION(fidhgr,dimID_y,len=my_GLO)
call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidhgr,dimID_x,len=mx_GLO)
call erreur(status,.TRUE.,"inq_dim_x")

!--         
ALLOCATE(  glamtGLO(mx_GLO,my_GLO) , gphitGLO(mx_GLO,my_GLO) , zglamtGLO(mx_GLO,my_GLO)  )
ALLOCATE(  glamfGLO(mx_GLO,my_GLO) , gphifGLO(mx_GLO,my_GLO) , zglamfGLO(mx_GLO,my_GLO)  )
ALLOCATE(  glamuGLO(mx_GLO,my_GLO) , gphiuGLO(mx_GLO,my_GLO) , zglamuGLO(mx_GLO,my_GLO)  )
ALLOCATE(  glamvGLO(mx_GLO,my_GLO) , gphivGLO(mx_GLO,my_GLO) , zglamvGLO(mx_GLO,my_GLO)  )

status = NF90_INQ_VARID(fidhgr,"glamt",glamt_ID) ; call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidhgr,"gphit",gphit_ID) ; call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidhgr,"glamf",glamf_ID) ; call erreur(status,.TRUE.,"inq_glamf_ID")
status = NF90_INQ_VARID(fidhgr,"gphif",gphif_ID) ; call erreur(status,.TRUE.,"inq_gphif_ID")
status = NF90_INQ_VARID(fidhgr,"glamu",glamu_ID) ; call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidhgr,"gphiu",gphiu_ID) ; call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidhgr,"glamv",glamv_ID) ; call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidhgr,"gphiv",gphiv_ID) ; call erreur(status,.TRUE.,"inq_gphiv_ID")

status = NF90_GET_VAR(fidhgr,glamt_ID,glamtGLO) ; call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidhgr,gphit_ID,gphitGLO) ; call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidhgr,glamf_ID,glamfGLO) ; call erreur(status,.TRUE.,"getvar_glamf")
status = NF90_GET_VAR(fidhgr,gphif_ID,gphifGLO) ; call erreur(status,.TRUE.,"getvar_gphif")
status = NF90_GET_VAR(fidhgr,glamu_ID,glamuGLO) ; call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidhgr,gphiu_ID,gphiuGLO) ; call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidhgr,glamv_ID,glamvGLO) ; call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidhgr,gphiv_ID,gphivGLO) ; call erreur(status,.TRUE.,"getvar_gphiv")                                          

status = NF90_CLOSE(fidhgr)                      
call erreur(status,.TRUE.,"end_read_horizontal_grid") 

!=========================================================================================
! 2- Read the grid coordinates for the REGIONAL domain
!=========================================================================================

write(*,*) 'Reading ', TRIM(file_in_coord_REG)
                                             
status = NF90_OPEN(TRIM(file_in_coord_REG),0,fidgridout)          
call erreur(status,.TRUE.,"read_output_grid") 
                                                     
 status = NF90_INQ_DIMID(fidgridout,"y",dimID_y)
 call erreur(status,.TRUE.,"inq_dimID_y")
 status = NF90_INQ_DIMID(fidgridout,"x",dimID_x)
 call erreur(status,.TRUE.,"inq_dimID_x")
                                                     
 status = NF90_INQUIRE_DIMENSION(fidgridout,dimID_y,len=my_REG)
 call erreur(status,.TRUE.,"inq_dim_y")
 status = NF90_INQUIRE_DIMENSION(fidgridout,dimID_x,len=mx_REG)
 call erreur(status,.TRUE.,"inq_dim_x")

 ALLOCATE(  gphitREG(mx_REG,my_REG)  )
 ALLOCATE(  glamtREG(mx_REG,my_REG)  )
 ALLOCATE(  zglamtREG(mx_REG,my_REG) )
 ALLOCATE(  gphifREG(mx_REG,my_REG)  )
 ALLOCATE(  glamfREG(mx_REG,my_REG)  )
 ALLOCATE(  zglamfREG(mx_REG,my_REG) )
 ALLOCATE(  gphiuREG(mx_REG,my_REG)  )
 ALLOCATE(  glamuREG(mx_REG,my_REG)  )
 ALLOCATE(  zglamuREG(mx_REG,my_REG) )
 ALLOCATE(  gphivREG(mx_REG,my_REG)  )
 ALLOCATE(  glamvREG(mx_REG,my_REG)  )
 ALLOCATE(  zglamvREG(mx_REG,my_REG) )
 ALLOCATE(  angoutXt(mx_REG,my_REG) , angleXt(mx_REG,my_REG) )
 ALLOCATE(  angoutYt(mx_REG,my_REG) , angleYt(mx_REG,my_REG) )
 ALLOCATE(  angoutXf(mx_REG,my_REG) , angleXf(mx_REG,my_REG) )
 ALLOCATE(  angoutYf(mx_REG,my_REG) , angleYf(mx_REG,my_REG) )
 ALLOCATE(  angoutXu(mx_REG,my_REG) , angleXu(mx_REG,my_REG) )
 ALLOCATE(  angoutYu(mx_REG,my_REG) , angleYu(mx_REG,my_REG) )
 ALLOCATE(  angoutXv(mx_REG,my_REG) , angleXv(mx_REG,my_REG) )
 ALLOCATE(  angoutYv(mx_REG,my_REG) , angleYv(mx_REG,my_REG) )

 status = NF90_INQ_VARID(fidgridout,"gphit",gphitREG_ID)     ; call erreur(status,.TRUE.,"inq_gphitREG_ID")
 status = NF90_INQ_VARID(fidgridout,"glamt",glamtREG_ID)     ; call erreur(status,.TRUE.,"inq_glamtREG_ID")
 status = NF90_INQ_VARID(fidgridout,"gphif",gphifREG_ID)     ; call erreur(status,.TRUE.,"inq_gphifREG_ID")
 status = NF90_INQ_VARID(fidgridout,"glamf",glamfREG_ID)     ; call erreur(status,.TRUE.,"inq_glamfREG_ID")
 status = NF90_INQ_VARID(fidgridout,"gphiu",gphiuREG_ID)     ; call erreur(status,.TRUE.,"inq_gphiuREG_ID")
 status = NF90_INQ_VARID(fidgridout,"glamu",glamuREG_ID)     ; call erreur(status,.TRUE.,"inq_glamuREG_ID")
 status = NF90_INQ_VARID(fidgridout,"gphiv",gphivREG_ID)     ; call erreur(status,.TRUE.,"inq_gphivREG_ID")
 status = NF90_INQ_VARID(fidgridout,"glamv",glamvREG_ID)     ; call erreur(status,.TRUE.,"inq_glamvREG_ID")
                                                    
 status = NF90_GET_VAR(fidgridout,gphitREG_ID,gphitREG)     ; call erreur(status,.TRUE.,"getvar_gphitREG")
 status = NF90_GET_VAR(fidgridout,glamtREG_ID,glamtREG)     ; call erreur(status,.TRUE.,"getvar_glamtREG")
 status = NF90_GET_VAR(fidgridout,gphifREG_ID,gphifREG)     ; call erreur(status,.TRUE.,"getvar_gphifREG")
 status = NF90_GET_VAR(fidgridout,glamfREG_ID,glamfREG)     ; call erreur(status,.TRUE.,"getvar_glamfREG")
 status = NF90_GET_VAR(fidgridout,gphiuREG_ID,gphiuREG)     ; call erreur(status,.TRUE.,"getvar_gphiuREG")
 status = NF90_GET_VAR(fidgridout,glamuREG_ID,glamuREG)     ; call erreur(status,.TRUE.,"getvar_glamuREG")
 status = NF90_GET_VAR(fidgridout,gphivREG_ID,gphivREG)     ; call erreur(status,.TRUE.,"getvar_gphivREG")
 status = NF90_GET_VAR(fidgridout,glamvREG_ID,glamvREG)     ; call erreur(status,.TRUE.,"getvar_glamvREG") 
                                           
status = NF90_CLOSE(fidgridout)                      
call erreur(status,.TRUE.,"fin_lecture_output_grid")     

!-- Rotation angles 
if ( ln_dateline ) then
  where ( glamtREG(:,:) .lt. 0.0 )
    zglamtREG(:,:) = 360.0 + glamtREG(:,:)
  elsewhere
    zglamtREG(:,:) = glamtREG(:,:)
  endwhere
  !-
  where ( glamfREG(:,:) .lt. 0.0 )
    zglamfREG(:,:) = 360.0 + glamfREG(:,:)
  elsewhere
    zglamfREG(:,:) = glamfREG(:,:)
  endwhere
  !-
  where ( glamuREG(:,:) .lt. 0.0 )
    zglamuREG(:,:) = 360.0 + glamuREG(:,:)
  elsewhere
    zglamuREG(:,:) = glamuREG(:,:)
  endwhere
  !-
  where ( glamvREG(:,:) .lt. 0.0 )
    zglamvREG(:,:) = 360.0 + glamvREG(:,:)
  elsewhere
    zglamvREG(:,:) = glamvREG(:,:)
  endwhere
else
  zglamtREG(:,:) = glamtREG(:,:)
  zglamfREG(:,:) = glamfREG(:,:)
  zglamuREG(:,:) = glamuREG(:,:)
  zglamvREG(:,:) = glamvREG(:,:)
endif

do iREG=1,mx_REG
do jREG=1,my_REG

  ! local angle between i-direction and the zonal direction and between j-direction and the meridional direction (should be similar)

  if     ( iREG .eq. 1      ) then
   angoutXt(iREG,jREG) = ATAN2( gphitREG(iREG+1,jREG) - gphitREG(iREG  ,jREG) , ( zglamtREG(iREG+1,jREG) - zglamtREG(iREG  ,jREG) ) * cos(rad*gphitREG(iREG,jREG)) )
  elseif ( iREG .eq. mx_REG ) then 
   angoutXt(iREG,jREG) = ATAN2( gphitREG(iREG  ,jREG) - gphitREG(iREG-1,jREG) , ( zglamtREG(iREG  ,jREG) - zglamtREG(iREG-1,jREG) ) * cos(rad*gphitREG(iREG,jREG)) )
  else
   angoutXt(iREG,jREG) = ATAN2( gphitREG(iREG+1,jREG) - gphitREG(iREG-1,jREG) , ( zglamtREG(iREG+1,jREG) - zglamtREG(iREG-1,jREG) ) * cos(rad*gphitREG(iREG,jREG)) )
  endif
  if     ( jREG .eq. 1      ) then
   angoutYt(iREG,jREG) = ATAN2( gphitREG(iREG,jREG+1) - gphitREG(iREG,jREG  ) , ( zglamtREG(iREG,jREG+1) - zglamtREG(iREG,jREG  ) ) * cos(rad*gphitREG(iREG,jREG)) )
  elseif ( jREG .eq. my_REG ) then
   angoutYt(iREG,jREG) = ATAN2( gphitREG(iREG,jREG  ) - gphitREG(iREG,jREG-1) , ( zglamtREG(iREG,jREG  ) - zglamtREG(iREG,jREG-1) ) * cos(rad*gphitREG(iREG,jREG)) )
  else
   angoutYt(iREG,jREG) = ATAN2( gphitREG(iREG,jREG+1) - gphitREG(iREG,jREG-1) , ( zglamtREG(iREG,jREG+1) - zglamtREG(iREG,jREG-1) ) * cos(rad*gphitREG(iREG,jREG)) )
  endif

  if     ( iREG .eq. 1      ) then
   angoutXf(iREG,jREG) = ATAN2( gphifREG(iREG+1,jREG) - gphifREG(iREG  ,jREG) , ( zglamfREG(iREG+1,jREG) - zglamfREG(iREG  ,jREG) ) * cos(rad*gphifREG(iREG,jREG)) )
  elseif ( iREG .eq. mx_REG ) then 
   angoutXf(iREG,jREG) = ATAN2( gphifREG(iREG  ,jREG) - gphifREG(iREG-1,jREG) , ( zglamfREG(iREG  ,jREG) - zglamfREG(iREG-1,jREG) ) * cos(rad*gphifREG(iREG,jREG)) )
  else
   angoutXf(iREG,jREG) = ATAN2( gphifREG(iREG+1,jREG) - gphifREG(iREG-1,jREG) , ( zglamfREG(iREG+1,jREG) - zglamfREG(iREG-1,jREG) ) * cos(rad*gphifREG(iREG,jREG)) )
  endif
  if     ( jREG .eq. 1      ) then
   angoutYf(iREG,jREG) = ATAN2( gphifREG(iREG,jREG+1) - gphifREG(iREG,jREG  ) , ( zglamfREG(iREG,jREG+1) - zglamfREG(iREG,jREG  ) ) * cos(rad*gphifREG(iREG,jREG)) )
  elseif ( jREG .eq. my_REG ) then
   angoutYf(iREG,jREG) = ATAN2( gphifREG(iREG,jREG  ) - gphifREG(iREG,jREG-1) , ( zglamfREG(iREG,jREG  ) - zglamfREG(iREG,jREG-1) ) * cos(rad*gphifREG(iREG,jREG)) )
  else
   angoutYf(iREG,jREG) = ATAN2( gphifREG(iREG,jREG+1) - gphifREG(iREG,jREG-1) , ( zglamfREG(iREG,jREG+1) - zglamfREG(iREG,jREG-1) ) * cos(rad*gphifREG(iREG,jREG)) )
  endif

  if     ( iREG .eq. 1      ) then
   angoutXu(iREG,jREG) = ATAN2( gphiuREG(iREG+1,jREG) - gphiuREG(iREG  ,jREG) , ( zglamuREG(iREG+1,jREG) - zglamuREG(iREG  ,jREG) ) * cos(rad*gphiuREG(iREG,jREG)) )
  elseif ( iREG .eq. mx_REG ) then 
   angoutXu(iREG,jREG) = ATAN2( gphiuREG(iREG  ,jREG) - gphiuREG(iREG-1,jREG) , ( zglamuREG(iREG  ,jREG) - zglamuREG(iREG-1,jREG) ) * cos(rad*gphiuREG(iREG,jREG)) )
  else
   angoutXu(iREG,jREG) = ATAN2( gphiuREG(iREG+1,jREG) - gphiuREG(iREG-1,jREG) , ( zglamuREG(iREG+1,jREG) - zglamuREG(iREG-1,jREG) ) * cos(rad*gphiuREG(iREG,jREG)) )
  endif
  if     ( jREG .eq. 1      ) then
   angoutYu(iREG,jREG) = ATAN2( gphiuREG(iREG,jREG+1) - gphiuREG(iREG,jREG  ) , ( zglamuREG(iREG,jREG+1) - zglamuREG(iREG,jREG  ) ) * cos(rad*gphiuREG(iREG,jREG)) )
  elseif ( jREG .eq. my_REG ) then
   angoutYu(iREG,jREG) = ATAN2( gphiuREG(iREG,jREG  ) - gphiuREG(iREG,jREG-1) , ( zglamuREG(iREG,jREG  ) - zglamuREG(iREG,jREG-1) ) * cos(rad*gphiuREG(iREG,jREG)) )
  else
   angoutYu(iREG,jREG) = ATAN2( gphiuREG(iREG,jREG+1) - gphiuREG(iREG,jREG-1) , ( zglamuREG(iREG,jREG+1) - zglamuREG(iREG,jREG-1) ) * cos(rad*gphiuREG(iREG,jREG)) )
  endif

  if     ( iREG .eq. 1      ) then
   angoutXv(iREG,jREG) = ATAN2( gphivREG(iREG+1,jREG) - gphivREG(iREG  ,jREG) , ( zglamvREG(iREG+1,jREG) - zglamvREG(iREG  ,jREG) ) * cos(rad*gphivREG(iREG,jREG)) )
  elseif ( iREG .eq. mx_REG ) then 
   angoutXv(iREG,jREG) = ATAN2( gphivREG(iREG  ,jREG) - gphivREG(iREG-1,jREG) , ( zglamvREG(iREG  ,jREG) - zglamvREG(iREG-1,jREG) ) * cos(rad*gphivREG(iREG,jREG)) )
  else
   angoutXv(iREG,jREG) = ATAN2( gphivREG(iREG+1,jREG) - gphivREG(iREG-1,jREG) , ( zglamvREG(iREG+1,jREG) - zglamvREG(iREG-1,jREG) ) * cos(rad*gphivREG(iREG,jREG)) )
  endif
  if     ( jREG .eq. 1      ) then
   angoutYv(iREG,jREG) = ATAN2( gphivREG(iREG,jREG+1) - gphivREG(iREG,jREG  ) , ( zglamvREG(iREG,jREG+1) - zglamvREG(iREG,jREG  ) ) * cos(rad*gphivREG(iREG,jREG)) )
  elseif ( jREG .eq. my_REG ) then
   angoutYv(iREG,jREG) = ATAN2( gphivREG(iREG,jREG  ) - gphivREG(iREG,jREG-1) , ( zglamvREG(iREG,jREG  ) - zglamvREG(iREG,jREG-1) ) * cos(rad*gphivREG(iREG,jREG)) )
  else
   angoutYv(iREG,jREG) = ATAN2( gphivREG(iREG,jREG+1) - gphivREG(iREG,jREG-1) , ( zglamvREG(iREG,jREG+1) - zglamvREG(iREG,jREG-1) ) * cos(rad*gphivREG(iREG,jREG)) )
  endif

enddo
enddo

!===============

!=========================================================================================
! 3- CALCULATE WEIGHTS TO INTERPOLATE GLOBAL FILEDS ONTO THE REGIONAL GRID 
!    (and find what are the 4 global T-points surrounding a T-point of the boundary)
!=========================================================================================

ALLOCATE( ziWt  (mx_REG,my_REG), wgNEt(mx_REG,my_REG) )
ALLOCATE( ziEt  (mx_REG,my_REG), wgSEt(mx_REG,my_REG) )
ALLOCATE( zjSt  (mx_REG,my_REG), wgSWt(mx_REG,my_REG) )
ALLOCATE( zjNt  (mx_REG,my_REG), wgNWt(mx_REG,my_REG) )
ALLOCATE( zzit  (mx_REG,my_REG), zzjt  (mx_REG,my_REG))
ALLOCATE( zzitm1(mx_REG,my_REG), zzjtm1(mx_REG,my_REG))
ALLOCATE( zzitp1(mx_REG,my_REG), zzjtp1(mx_REG,my_REG))
!-
ALLOCATE( ziWf  (mx_REG,my_REG), wgNEf(mx_REG,my_REG) )
ALLOCATE( ziEf  (mx_REG,my_REG), wgSEf(mx_REG,my_REG) )
ALLOCATE( zjSf  (mx_REG,my_REG), wgSWf(mx_REG,my_REG) )
ALLOCATE( zjNf  (mx_REG,my_REG), wgNWf(mx_REG,my_REG) )
ALLOCATE( zzif  (mx_REG,my_REG), zzjf  (mx_REG,my_REG))
ALLOCATE( zzifm1(mx_REG,my_REG), zzjfm1(mx_REG,my_REG))
ALLOCATE( zzifp1(mx_REG,my_REG), zzjfp1(mx_REG,my_REG))
!-
ALLOCATE( ziWu  (mx_REG,my_REG), wgNEu(mx_REG,my_REG) )
ALLOCATE( ziEu  (mx_REG,my_REG), wgSEu(mx_REG,my_REG) )
ALLOCATE( zjSu  (mx_REG,my_REG), wgSWu(mx_REG,my_REG) )
ALLOCATE( zjNu  (mx_REG,my_REG), wgNWu(mx_REG,my_REG) )
ALLOCATE( zziu  (mx_REG,my_REG), zzju  (mx_REG,my_REG))
ALLOCATE( zzium1(mx_REG,my_REG), zzjum1(mx_REG,my_REG))
ALLOCATE( zziup1(mx_REG,my_REG), zzjup1(mx_REG,my_REG))
!-
ALLOCATE( ziWv  (mx_REG,my_REG), wgNEv(mx_REG,my_REG) )
ALLOCATE( ziEv  (mx_REG,my_REG), wgSEv(mx_REG,my_REG) )
ALLOCATE( zjSv  (mx_REG,my_REG), wgSWv(mx_REG,my_REG) )
ALLOCATE( zjNv  (mx_REG,my_REG), wgNWv(mx_REG,my_REG) )
ALLOCATE( zziv  (mx_REG,my_REG), zzjv  (mx_REG,my_REG))
ALLOCATE( zzivm1(mx_REG,my_REG), zzjvm1(mx_REG,my_REG))
ALLOCATE( zzivp1(mx_REG,my_REG), zzjvp1(mx_REG,my_REG))

IF ( ln_dateline ) THEN
  where ( glamtGLO(:,:) .lt. 0.0 )
    zglamtGLO(:,:) = 360.0 + glamtGLO(:,:)
  elsewhere
    zglamtGLO(:,:) = glamtGLO(:,:)
  endwhere
  !-
  where ( glamfGLO(:,:) .lt. 0.0 )
    zglamfGLO(:,:) = 360.0 + glamfGLO(:,:)
  elsewhere
    zglamfGLO(:,:) = glamfGLO(:,:)
  endwhere
  !-
  where ( glamuGLO(:,:) .lt. 0.0 )
    zglamuGLO(:,:) = 360.0 + glamuGLO(:,:)
  elsewhere
    zglamuGLO(:,:) = glamuGLO(:,:)
  endwhere
  !-
  where ( glamvGLO(:,:) .lt. 0.0 )
    zglamvGLO(:,:) = 360.0 + glamvGLO(:,:)
  elsewhere
    zglamvGLO(:,:) = glamvGLO(:,:)
  endwhere
ELSE
  zglamtGLO(:,:) = glamtGLO(:,:) 
  zglamfGLO(:,:) = glamfGLO(:,:) 
  zglamuGLO(:,:) = glamuGLO(:,:)
  zglamvGLO(:,:) = glamvGLO(:,:)
ENDIF

zzit(:,:)=0 !- to notice error if a point is not filled
zzjt(:,:)=0 !- to notice error if a point is not filled

!##### gridT #####
write(*,*) 'Coefficients for T points...'

DO iREG=1,mx_REG
DO jREG=1,my_REG

  distmin = ra * rad * 180.0

  !-- loop to find closest point
  do iGLO=1,mx_GLO
  do jGLO=1,my_GLO
    dist = ra * rad * sqrt(   ( cos(rad*gphitREG(iREG,jREG)) * ( zglamtREG(iREG,jREG) - zglamtGLO(iGLO,jGLO) ) )**2  &
      &                     + (                                    gphitREG(iREG,jREG) -  gphitGLO(iGLO,jGLO)   )**2  )
    if ( dist .lt. distmin ) then
      distmin=dist
      zzit(iREG,jREG) = iGLO
      zzjt(iREG,jREG) = jGLO
    endif
  enddo 
  enddo
  !- East-West bounds
  if ( nn_perio .eq. 1 ) then !- periodic
    if     ( zzit(iREG,jREG)+1 .gt. mx_GLO ) then
      zzitp1(iREG,jREG) = 3
      zzitm1(iREG,jREG) = zzit(iREG,jREG) - 1
    elseif ( zzit(iREG,jREG)-1 .lt. 1     ) then
      zzitp1(iREG,jREG) = zzit(iREG,jREG) + 1
      zzitm1(iREG,jREG) = mx_GLO-2
    else
      zzitp1(iREG,jREG) = zzit(iREG,jREG) + 1
      zzitm1(iREG,jREG) = zzit(iREG,jREG) - 1
    endif
  elseif ( nn_perio .eq. 0 ) then
    if     ( zzit(iREG,jREG)+1 .gt. mx_GLO ) then
      zzitp1(iREG,jREG) = zzit(iREG,jREG)
      zzitm1(iREG,jREG) = zzit(iREG,jREG) - 1
    elseif ( zzit(iREG,jREG)-1 .lt. 1     ) then
      zzitp1(iREG,jREG) = zzit(iREG,jREG) + 1
      zzitm1(iREG,jREG) = zzit(iREG,jREG)
    else
      zzitp1(iREG,jREG) = zzit(iREG,jREG) + 1
      zzitm1(iREG,jREG) = zzit(iREG,jREG) - 1
    endif
  else
    write(*,*) '~!@#$%^* nn_perio must be either 0 or 1 >>>>> stop !!'
    stop
  endif
  !- upper lower bounds of global input grid
  if     ( zzjt(iREG,jREG)+1 .gt. my_GLO ) then
    zzjtp1(iREG,jREG) = zzjt(iREG,jREG)
    zzjtm1(iREG,jREG) = zzjt(iREG,jREG) - 1
  elseif ( zzit(iREG,jREG)-1 .lt. 1     ) then
    zzjtp1(iREG,jREG) = zzjt(iREG,jREG) + 1
    zzjtm1(iREG,jREG) = zzjt(iREG,jREG)
  else
    zzjtp1(iREG,jREG) = zzjt(iREG,jREG) + 1
    zzjtm1(iREG,jREG) = zzjt(iREG,jREG) - 1
  endif
  !-
  zzitm1(iREG,jREG) = MAX( MIN( zzitm1(iREG,jREG), mx_GLO ), 1 )
  zzitp1(iREG,jREG) = MAX( MIN( zzitp1(iREG,jREG), mx_GLO ), 1 )
  zzjtm1(iREG,jREG) = MAX( MIN( zzjtm1(iREG,jREG), my_GLO ), 1 )
  zzjtp1(iREG,jREG) = MAX( MIN( zzjtp1(iREG,jREG), my_GLO ), 1 )

  !-- local angle on global input grid between i-direction and the zonal direction and between j-direction and the meridional direction
  anginX = ATAN2(    gphitGLO( zzitp1(iREG,jREG) , zzjt(iREG,jREG) ) -  gphitGLO( zzitm1(iREG,jREG) , zzjt(iREG,jREG)) , &
  &                ( zglamtGLO( zzitp1(iREG,jREG) , zzjt(iREG,jREG) ) - zglamtGLO( zzitm1(iREG,jREG) , zzjt(iREG,jREG)) ) &
  &                * cos( gphitGLO( zzit(iREG,jREG) , zzjt(iREG,jREG) ) * rad )                                                )

  anginY = ATAN2(    gphitGLO( zzit(iREG,jREG) , zzjtp1(iREG,jREG) ) -  gphitGLO( zzit(iREG,jREG) , zzjtm1(iREG,jREG)) , &
  &                ( zglamtGLO( zzit(iREG,jREG) , zzjtp1(iREG,jREG) ) - zglamtGLO( zzit(iREG,jREG) , zzjtm1(iREG,jREG)) ) &
  &                * cos( gphitGLO( zzit(iREG,jREG) , zzjt(iREG,jREG) ) * rad )                                                )

  !-- local angle between the two grids :
  angleXt(iREG,jREG) = angoutXt(iREG,jREG) - anginX
  angleYt(iREG,jREG) = angoutYt(iREG,jREG) - anginY
  if ( abs(angleXt(iREG,jREG)/rad) .gt. 45.0 .or. abs(angleYt(iREG,jREG)/rad) .gt. 45.0 ) then
    write(*,*) '@@@@@@ WARNING : angle between 2 grids > 45 deg at point of the bdy :', iREG, jREG
  endif     

  !------------------------------------------------------------
  if ( zglamtREG(iREG,jREG) .ge. zglamtGLO(zzit(iREG,jREG),zzjt(iREG,jREG)) )  then
    ziWt(iREG,jREG)   = zzit  (iREG,jREG) ! closest point westward
    ziEt(iREG,jREG)   = zzitp1(iREG,jREG) ! closest point eastward
  else
    ziWt(iREG,jREG)   = zzitm1(iREG,jREG) ! closest point westward
    ziEt(iREG,jREG)   = zzit  (iREG,jREG) ! closest point eastward
  endif
  !--
  if ( gphitREG(iREG,jREG) .ge. gphitGLO(zzit(iREG,jREG),zzjt(iREG,jREG)) )  then
    zjSt(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point southward
    zjNt(iREG,jREG)   = zzjtp1(iREG,jREG) ! closest point northward
  else
    zjSt(iREG,jREG)   = zzjtm1(iREG,jREG) ! closest point southward
    zjNt(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point northward
  endif
  !--
  distSW = sqrt(   ( cos(rad*gphitREG(iREG,jREG)) * ( zglamtREG(iREG,jREG) - zglamtGLO(ziWt(iREG,jREG),zjSt(iREG,jREG)) ) )**2  &
      &          + (                                     gphitREG(iREG,jREG) -  gphitGLO(ziWt(iREG,jREG),zjSt(iREG,jREG)) )**2  )
  distSE = sqrt(   ( cos(rad*gphitREG(iREG,jREG)) * ( zglamtREG(iREG,jREG) - zglamtGLO(ziEt(iREG,jREG),zjSt(iREG,jREG)) ) )**2  &
      &          + (                                     gphitREG(iREG,jREG) -  gphitGLO(ziEt(iREG,jREG),zjSt(iREG,jREG)) )**2  )
  distNW = sqrt(   ( cos(rad*gphitREG(iREG,jREG)) * ( zglamtREG(iREG,jREG) - zglamtGLO(ziWt(iREG,jREG),zjNt(iREG,jREG)) ) )**2  &
      &          + (                                     gphitREG(iREG,jREG) -  gphitGLO(ziWt(iREG,jREG),zjNt(iREG,jREG)) )**2  )
  distNE = sqrt(   ( cos(rad*gphitREG(iREG,jREG)) * ( zglamtREG(iREG,jREG) - zglamtGLO(ziEt(iREG,jREG),zjNt(iREG,jREG)) ) )**2  &
      &          + (                                     gphitREG(iREG,jREG) -  gphitGLO(ziEt(iREG,jREG),zjNt(iREG,jREG)) )**2  )
  !--
  wgSWt(iREG,jREG) =   ( zglamtGLO(ziEt(iREG,jREG),zjNt(iREG,jREG)) - zglamtREG(iREG,jREG) ) &
  &                  * (  gphitGLO(ziEt(iREG,jREG),zjNt(iREG,jREG)) -  gphitREG(iREG,jREG) )
  wgSEt(iREG,jREG) =   ( zglamtGLO(ziWt(iREG,jREG),zjNt(iREG,jREG)) - zglamtREG(iREG,jREG) ) &
  &                  * (  gphitREG(iREG,jREG) -  gphitGLO(ziWt(iREG,jREG),zjNt(iREG,jREG)) )
  wgNWt(iREG,jREG) =   ( zglamtREG(iREG,jREG) - zglamtGLO(ziEt(iREG,jREG),zjSt(iREG,jREG)) ) &
  &                  * (  gphitGLO(ziEt(iREG,jREG),zjSt(iREG,jREG)) -  gphitREG(iREG,jREG) )
  wgNEt(iREG,jREG) =   ( zglamtREG(iREG,jREG) - zglamtGLO(ziWt(iREG,jREG),zjSt(iREG,jREG)) ) &
  &                  * (  gphitREG(iREG,jREG) -  gphitGLO(ziWt(iREG,jREG),zjSt(iREG,jREG)) )

ENDDO !-- jREG
ENDDO !-- iREG

!##### gridF #####
write(*,*) 'Coefficients for F points...'

zzif(:,:)=0 !- to notice error if a point is not filled
zzjf(:,:)=0 !- to notice error if a point is not filled

DO iREG=1,mx_REG
DO jREG=1,my_REG

  distmin = ra * rad * 180.0

  !-- loop to find closest point
  do iGLO=1,mx_GLO
  do jGLO=1,my_GLO
    dist = ra * rad * sqrt(   ( cos(rad*gphifREG(iREG,jREG)) * ( zglamfREG(iREG,jREG) - zglamfGLO(iGLO,jGLO) ) )**2  &
      &                     + (                                   gphifREG(iREG,jREG) -  gphifGLO(iGLO,jGLO)   )**2  )
    if ( dist .lt. distmin ) then
      distmin=dist
      zzif(iREG,jREG) = iGLO
      zzjf(iREG,jREG) = jGLO
    endif
  enddo 
  enddo
  !- East-West bounds
  if ( nn_perio .eq. 1 ) then !- periodic
    if     ( zzif(iREG,jREG)+1 .gt. mx_GLO ) then
      zzifp1(iREG,jREG) = 3 
      zzifm1(iREG,jREG) = zzif(iREG,jREG) - 1
    elseif ( zzif(iREG,jREG)-1 .lt. 1     ) then
      zzifp1(iREG,jREG) = zzif(iREG,jREG) + 1
      zzifm1(iREG,jREG) = mx_GLO-2 
    else
      zzifp1(iREG,jREG) = zzif(iREG,jREG) + 1
      zzifm1(iREG,jREG) = zzif(iREG,jREG) - 1
    endif
  else !- non-periodic
    if     ( zzif(iREG,jREG)+1 .gt. mx_GLO ) then
      zzifp1(iREG,jREG) = zzif(iREG,jREG)
      zzifm1(iREG,jREG) = zzif(iREG,jREG) - 1
    elseif ( zzif(iREG,jREG)-1 .lt. 1     ) then
      zzifp1(iREG,jREG) = zzif(iREG,jREG) + 1
      zzifm1(iREG,jREG) = zzif(iREG,jREG)
    else
      zzifp1(iREG,jREG) = zzif(iREG,jREG) + 1
      zzifm1(iREG,jREG) = zzif(iREG,jREG) - 1
    endif
  endif
  !- upper lower bounds of global input grid
  if     ( zzjf(iREG,jREG)+1 .gt. my_GLO ) then
    zzjfp1(iREG,jREG) = zzjf(iREG,jREG)
    zzjfm1(iREG,jREG) = zzjf(iREG,jREG) - 1
  elseif ( zzif(iREG,jREG)-1 .lt. 1     ) then
    zzjfp1(iREG,jREG) = zzjf(iREG,jREG) + 1
    zzjfm1(iREG,jREG) = zzjf(iREG,jREG)
  else
    zzjfp1(iREG,jREG) = zzjf(iREG,jREG) + 1
    zzjfm1(iREG,jREG) = zzjf(iREG,jREG) - 1
  endif
  !-
  zzifm1(iREG,jREG) = MAX( MIN( zzifm1(iREG,jREG), mx_GLO ), 1 )
  zzifp1(iREG,jREG) = MAX( MIN( zzifp1(iREG,jREG), mx_GLO ), 1 )
  zzjfm1(iREG,jREG) = MAX( MIN( zzjfm1(iREG,jREG), my_GLO ), 1 )
  zzjfp1(iREG,jREG) = MAX( MIN( zzjfp1(iREG,jREG), my_GLO ), 1 )

  !-- local angle on global input grid between i-direction and the zonal direction and between j-direction and the meridional direction
  anginX = ATAN2(    gphifGLO( zzifp1(iREG,jREG) , zzjf(iREG,jREG) ) -  gphifGLO( zzifm1(iREG,jREG) , zzjf(iREG,jREG)) , &
  &                ( zglamfGLO( zzifp1(iREG,jREG) , zzjf(iREG,jREG) ) - zglamfGLO( zzifm1(iREG,jREG) , zzjf(iREG,jREG)) ) &
  &                * cos( gphifGLO( zzif(iREG,jREG) , zzjf(iREG,jREG) ) * rad )                                                )

  anginY = ATAN2(    gphifGLO( zzif(iREG,jREG) , zzjfp1(iREG,jREG) ) -  gphifGLO( zzif(iREG,jREG) , zzjfm1(iREG,jREG)) , &
  &                ( zglamfGLO( zzif(iREG,jREG) , zzjfp1(iREG,jREG) ) - zglamfGLO( zzif(iREG,jREG) , zzjfm1(iREG,jREG)) ) &
  &                * cos( gphifGLO( zzif(iREG,jREG) , zzjf(iREG,jREG) ) * rad )                                                )

  !-- local angle between the two grids :
  angleXf(iREG,jREG) = angoutXf(iREG,jREG) - anginX
  angleYf(iREG,jREG) = angoutYt(iREG,jREG) - anginY
  if ( abs(angleXf(iREG,jREG)/rad) .gt. 45.0 .or. abs(angleYf(iREG,jREG)/rad) .gt. 45.0 ) then
    write(*,*) '@@@@@@ WARNING : angle between 2 grids > 45 deg at point of the bdy :', iREG, jREG
  endif     

  !------------------------------------------------------------
  if ( zglamfREG(iREG,jREG) .ge. zglamfGLO(zzif(iREG,jREG),zzjt(iREG,jREG)) )  then
    ziWf(iREG,jREG)   = zzif  (iREG,jREG) ! closest point westward
    ziEf(iREG,jREG)   = zzifp1(iREG,jREG) ! closest point eastward
  else
    ziWf(iREG,jREG)   = zzifm1(iREG,jREG) ! closest point westward
    ziEf(iREG,jREG)   = zzif  (iREG,jREG) ! closest point eastward
  endif
  !--
  if ( gphifREG(iREG,jREG) .ge. gphifGLO(zzif(iREG,jREG),zzjt(iREG,jREG)) )  then
    zjSf(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point southward
    zjNf(iREG,jREG)   = zzjtp1(iREG,jREG) ! closest point northward
  else
    zjSf(iREG,jREG)   = zzjtm1(iREG,jREG) ! closest point southward
    zjNf(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point northward
  endif
  !--
  distSW = sqrt(   ( cos(rad*gphifREG(iREG,jREG)) * ( zglamfREG(iREG,jREG) - zglamfGLO(ziWf(iREG,jREG),zjSf(iREG,jREG)) ) )**2  &
      &          + (                                     gphifREG(iREG,jREG) -  gphifGLO(ziWf(iREG,jREG),zjSf(iREG,jREG)) )**2  )
  distSE = sqrt(   ( cos(rad*gphifREG(iREG,jREG)) * ( zglamfREG(iREG,jREG) - zglamfGLO(ziEf(iREG,jREG),zjSf(iREG,jREG)) ) )**2  &
      &          + (                                     gphifREG(iREG,jREG) -  gphifGLO(ziEf(iREG,jREG),zjSf(iREG,jREG)) )**2  )
  distNW = sqrt(   ( cos(rad*gphifREG(iREG,jREG)) * ( zglamfREG(iREG,jREG) - zglamfGLO(ziWf(iREG,jREG),zjNf(iREG,jREG)) ) )**2  &
      &          + (                                     gphifREG(iREG,jREG) -  gphifGLO(ziWf(iREG,jREG),zjNf(iREG,jREG)) )**2  )
  distNE = sqrt(   ( cos(rad*gphifREG(iREG,jREG)) * ( zglamfREG(iREG,jREG) - zglamfGLO(ziEf(iREG,jREG),zjNf(iREG,jREG)) ) )**2  &
      &          + (                                     gphifREG(iREG,jREG) -  gphifGLO(ziEf(iREG,jREG),zjNf(iREG,jREG)) )**2  )
  !--
  wgSWf(iREG,jREG) =   ( zglamfGLO(ziEf(iREG,jREG),zjNf(iREG,jREG)) - zglamfREG(iREG,jREG) ) &
  &                  * (  gphifGLO(ziEf(iREG,jREG),zjNf(iREG,jREG)) -  gphifREG(iREG,jREG) )
  wgSEf(iREG,jREG) =   ( zglamfGLO(ziWf(iREG,jREG),zjNf(iREG,jREG)) - zglamfREG(iREG,jREG) ) &
  &                  * (  gphifREG(iREG,jREG) -  gphifGLO(ziWf(iREG,jREG),zjNf(iREG,jREG)) )
  wgNWf(iREG,jREG) =   ( zglamfREG(iREG,jREG) - zglamfGLO(ziEf(iREG,jREG),zjSf(iREG,jREG)) ) &
  &                  * (  gphifGLO(ziEf(iREG,jREG),zjSf(iREG,jREG)) -  gphifREG(iREG,jREG) )
  wgNEf(iREG,jREG) =   ( zglamfREG(iREG,jREG) - zglamfGLO(ziWf(iREG,jREG),zjSf(iREG,jREG)) ) &
  &                  * (  gphifREG(iREG,jREG) -  gphifGLO(ziWf(iREG,jREG),zjSf(iREG,jREG)) )

ENDDO !-- jREG
ENDDO !-- iREG

!##### gridU #####
write(*,*) 'Coefficients for U points...'

zziu(:,:)=0 !- to notice error if a point is not filled
zzju(:,:)=0 !- to notice error if a point is not filled

DO iREG=1,mx_REG
DO jREG=1,my_REG

  distmin = ra * rad * 180.0

  !-- loop to find closest point
  do iGLO=1,mx_GLO
  do jGLO=1,my_GLO
    dist = ra * rad * sqrt(   ( cos(rad*gphiuREG(iREG,jREG)) * ( zglamuREG(iREG,jREG) - zglamuGLO(iGLO,jGLO) ) )**2  &
      &                     + (                                   gphiuREG(iREG,jREG) -  gphiuGLO(iGLO,jGLO)   )**2  )
    if ( dist .lt. distmin ) then
      distmin=dist
      zziu(iREG,jREG) = iGLO
      zzju(iREG,jREG) = jGLO
    endif
  enddo 
  enddo
  !- East-West bounds
  if ( nn_perio .eq. 1 ) then !- periodic
    if     ( zziu(iREG,jREG)+1 .gt. mx_GLO ) then
      zziup1(iREG,jREG) = 3
      zzium1(iREG,jREG) = zziu(iREG,jREG) - 1
    elseif ( zziu(iREG,jREG)-1 .lt. 1     ) then
      zziup1(iREG,jREG) = zziu(iREG,jREG) + 1
      zzium1(iREG,jREG) = mx_GLO-2
    else
      zziup1(iREG,jREG) = zziu(iREG,jREG) + 1
      zzium1(iREG,jREG) = zziu(iREG,jREG) - 1
    endif
  else !- non-periodic
    if     ( zziu(iREG,jREG)+1 .gt. mx_GLO ) then
      zziup1(iREG,jREG) = zziu(iREG,jREG)
      zzium1(iREG,jREG) = zziu(iREG,jREG) - 1
    elseif ( zziu(iREG,jREG)-1 .lt. 1     ) then
      zziup1(iREG,jREG) = zziu(iREG,jREG) + 1
      zzium1(iREG,jREG) = zziu(iREG,jREG)
    else
      zziup1(iREG,jREG) = zziu(iREG,jREG) + 1
      zzium1(iREG,jREG) = zziu(iREG,jREG) - 1
    endif
  endif
  !- upper lower bounds of global input grid
  if     ( zzju(iREG,jREG)+1 .gt. my_GLO ) then
    zzjup1(iREG,jREG) = zzju(iREG,jREG)
    zzjum1(iREG,jREG) = zzju(iREG,jREG) - 1
  elseif ( zziu(iREG,jREG)-1 .lt. 1     ) then
    zzjup1(iREG,jREG) = zzju(iREG,jREG) + 1
    zzjum1(iREG,jREG) = zzju(iREG,jREG)
  else
    zzjup1(iREG,jREG) = zzju(iREG,jREG) + 1
    zzjum1(iREG,jREG) = zzju(iREG,jREG) - 1
  endif
  !-
  zzium1(iREG,jREG) = MAX( MIN( zzium1(iREG,jREG), mx_GLO ), 1 )
  zziup1(iREG,jREG) = MAX( MIN( zziup1(iREG,jREG), mx_GLO ), 1 )
  zzjum1(iREG,jREG) = MAX( MIN( zzjum1(iREG,jREG), my_GLO ), 1 )
  zzjup1(iREG,jREG) = MAX( MIN( zzjup1(iREG,jREG), my_GLO ), 1 )

  !-- local angle on global input grid between i-direction and the zonal direction and between j-direction and the meridional direction
  anginX = ATAN2(    gphiuGLO( zziup1(iREG,jREG) , zzju(iREG,jREG) ) -  gphiuGLO( zzium1(iREG,jREG) , zzju(iREG,jREG)) , &
  &                ( zglamuGLO( zziup1(iREG,jREG) , zzju(iREG,jREG) ) - zglamuGLO( zzium1(iREG,jREG) , zzju(iREG,jREG)) ) &
  &                * cos( gphiuGLO( zziu(iREG,jREG) , zzju(iREG,jREG) ) * rad )                                                )

  anginY = ATAN2(    gphiuGLO( zziu(iREG,jREG) , zzjup1(iREG,jREG) ) -  gphiuGLO( zziu(iREG,jREG) , zzjum1(iREG,jREG)) , &
  &                ( zglamuGLO( zziu(iREG,jREG) , zzjup1(iREG,jREG) ) - zglamuGLO( zziu(iREG,jREG) , zzjum1(iREG,jREG)) ) &
  &                * cos( gphiuGLO( zziu(iREG,jREG) , zzju(iREG,jREG) ) * rad )                                                )

  !-- local angle between the two grids :
  angleXu(iREG,jREG) = angoutXu(iREG,jREG) - anginX
  angleYu(iREG,jREG) = angoutYt(iREG,jREG) - anginY
  if ( abs(angleXu(iREG,jREG)/rad) .gt. 45.0 .or. abs(angleYu(iREG,jREG)/rad) .gt. 45.0 ) then
    write(*,*) '@@@@@@ WARNING : angle between 2 grids > 45 deg at point of the bdy :', iREG, jREG
  endif     

  !------------------------------------------------------------
  if ( zglamuREG(iREG,jREG) .ge. zglamuGLO(zziu(iREG,jREG),zzjt(iREG,jREG)) )  then
    ziWu(iREG,jREG)   = zziu  (iREG,jREG) ! closest point westward
    ziEu(iREG,jREG)   = zziup1(iREG,jREG) ! closest point eastward
  else
    ziWu(iREG,jREG)   = zzium1(iREG,jREG) ! closest point westward
    ziEu(iREG,jREG)   = zziu  (iREG,jREG) ! closest point eastward
  endif
  !--
  if ( gphiuREG(iREG,jREG) .ge. gphiuGLO(zziu(iREG,jREG),zzjt(iREG,jREG)) )  then
    zjSu(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point southward
    zjNu(iREG,jREG)   = zzjtp1(iREG,jREG) ! closest point northward
  else
    zjSu(iREG,jREG)   = zzjtm1(iREG,jREG) ! closest point southward
    zjNu(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point northward
  endif
  !--
  distSW = sqrt(   ( cos(rad*gphiuREG(iREG,jREG)) * ( zglamuREG(iREG,jREG) - zglamuGLO(ziWu(iREG,jREG),zjSu(iREG,jREG)) ) )**2  &
      &          + (                                     gphiuREG(iREG,jREG) -  gphiuGLO(ziWu(iREG,jREG),zjSu(iREG,jREG)) )**2  )
  distSE = sqrt(   ( cos(rad*gphiuREG(iREG,jREG)) * ( zglamuREG(iREG,jREG) - zglamuGLO(ziEu(iREG,jREG),zjSu(iREG,jREG)) ) )**2  &
      &          + (                                     gphiuREG(iREG,jREG) -  gphiuGLO(ziEu(iREG,jREG),zjSu(iREG,jREG)) )**2  )
  distNW = sqrt(   ( cos(rad*gphiuREG(iREG,jREG)) * ( zglamuREG(iREG,jREG) - zglamuGLO(ziWu(iREG,jREG),zjNu(iREG,jREG)) ) )**2  &
      &          + (                                     gphiuREG(iREG,jREG) -  gphiuGLO(ziWu(iREG,jREG),zjNu(iREG,jREG)) )**2  )
  distNE = sqrt(   ( cos(rad*gphiuREG(iREG,jREG)) * ( zglamuREG(iREG,jREG) - zglamuGLO(ziEu(iREG,jREG),zjNu(iREG,jREG)) ) )**2  &
      &          + (                                     gphiuREG(iREG,jREG) -  gphiuGLO(ziEu(iREG,jREG),zjNu(iREG,jREG)) )**2  )
  !--
  wgSWu(iREG,jREG) =   ( zglamuGLO(ziEu(iREG,jREG),zjNu(iREG,jREG)) - zglamuREG(iREG,jREG) ) &
  &                  * (  gphiuGLO(ziEu(iREG,jREG),zjNu(iREG,jREG)) -  gphiuREG(iREG,jREG) )
  wgSEu(iREG,jREG) =   ( zglamuGLO(ziWu(iREG,jREG),zjNu(iREG,jREG)) - zglamuREG(iREG,jREG) ) &
  &                  * (  gphiuREG(iREG,jREG) -  gphiuGLO(ziWu(iREG,jREG),zjNu(iREG,jREG)) )
  wgNWu(iREG,jREG) =   ( zglamuREG(iREG,jREG) - zglamuGLO(ziEu(iREG,jREG),zjSu(iREG,jREG)) ) &
  &                  * (  gphiuGLO(ziEu(iREG,jREG),zjSu(iREG,jREG)) -  gphiuREG(iREG,jREG) )
  wgNEu(iREG,jREG) =   ( zglamuREG(iREG,jREG) - zglamuGLO(ziWu(iREG,jREG),zjSu(iREG,jREG)) ) &
  &                  * (  gphiuREG(iREG,jREG) -  gphiuGLO(ziWu(iREG,jREG),zjSu(iREG,jREG)) )

ENDDO !-- jREG
ENDDO !-- iREG

!##### gridV #####
write(*,*) 'Coefficients for V points...'

zziv(:,:)=0 !- to notice error if a point is not filled
zzjv(:,:)=0 !- to notice error if a point is not filled

DO iREG=1,mx_REG
DO jREG=1,my_REG

  distmin = ra * rad * 180.0

  !-- loop to find closest point
  do iGLO=1,mx_GLO
  do jGLO=1,my_GLO
    dist = ra * rad * sqrt(   ( cos(rad*gphivREG(iREG,jREG)) * ( zglamvREG(iREG,jREG) - zglamvGLO(iGLO,jGLO) ) )**2  &
      &                     + (                                   gphivREG(iREG,jREG) -  gphivGLO(iGLO,jGLO)   )**2  )
    if ( dist .lt. distmin ) then
      distmin=dist
      zziv(iREG,jREG) = iGLO
      zzjv(iREG,jREG) = jGLO
    endif
  enddo 
  enddo
  !- East-West bounds
  if ( nn_perio .eq. 1 ) then !- periodic
    if     ( zziv(iREG,jREG)+1 .gt. mx_GLO ) then
      zzivp1(iREG,jREG) = 3 
      zzivm1(iREG,jREG) = zziv(iREG,jREG) - 1
    elseif ( zziv(iREG,jREG)-1 .lt. 1     ) then
      zzivp1(iREG,jREG) = zziv(iREG,jREG) + 1
      zzivm1(iREG,jREG) = mx_GLO-2 
    else
      zzivp1(iREG,jREG) = zziv(iREG,jREG) + 1
      zzivm1(iREG,jREG) = zziv(iREG,jREG) - 1
    endif
  else !- non-periodic
    if     ( zziv(iREG,jREG)+1 .gt. mx_GLO ) then
      zzivp1(iREG,jREG) = zziv(iREG,jREG)
      zzivm1(iREG,jREG) = zziv(iREG,jREG) - 1
    elseif ( zziv(iREG,jREG)-1 .lt. 1     ) then
      zzivp1(iREG,jREG) = zziv(iREG,jREG) + 1
      zzivm1(iREG,jREG) = zziv(iREG,jREG)
    else
      zzivp1(iREG,jREG) = zziv(iREG,jREG) + 1
      zzivm1(iREG,jREG) = zziv(iREG,jREG) - 1
    endif
  endif
  !- upper lower bounds of global input grid
  if     ( zzjv(iREG,jREG)+1 .gt. my_GLO ) then
    zzjvp1(iREG,jREG) = zzjv(iREG,jREG)
    zzjvm1(iREG,jREG) = zzjv(iREG,jREG) - 1
  elseif ( zziv(iREG,jREG)-1 .lt. 1     ) then
    zzjvp1(iREG,jREG) = zzjv(iREG,jREG) + 1
    zzjvm1(iREG,jREG) = zzjv(iREG,jREG)
  else
    zzjvp1(iREG,jREG) = zzjv(iREG,jREG) + 1
    zzjvm1(iREG,jREG) = zzjv(iREG,jREG) - 1
  endif
  !-
  zzivm1(iREG,jREG) = MAX( MIN( zzivm1(iREG,jREG), mx_GLO ), 1 )
  zzivp1(iREG,jREG) = MAX( MIN( zzivp1(iREG,jREG), mx_GLO ), 1 )
  zzjvm1(iREG,jREG) = MAX( MIN( zzjvm1(iREG,jREG), my_GLO ), 1 )
  zzjvp1(iREG,jREG) = MAX( MIN( zzjvp1(iREG,jREG), my_GLO ), 1 )

  !-- local angle on global input grid between i-direction and the zonal direction and between j-direction and the meridional direction
  anginX = ATAN2(    gphivGLO( zzivp1(iREG,jREG) , zzjv(iREG,jREG) ) -  gphivGLO( zzivm1(iREG,jREG) , zzjv(iREG,jREG)) , &
  &                ( zglamvGLO( zzivp1(iREG,jREG) , zzjv(iREG,jREG) ) - zglamvGLO( zzivm1(iREG,jREG) , zzjv(iREG,jREG)) ) &
  &                * cos( gphivGLO( zziv(iREG,jREG) , zzjv(iREG,jREG) ) * rad )                                                )

  anginY = ATAN2(    gphivGLO( zziv(iREG,jREG) , zzjvp1(iREG,jREG) ) -  gphivGLO( zziv(iREG,jREG) , zzjvm1(iREG,jREG)) , &
  &                ( zglamvGLO( zziv(iREG,jREG) , zzjvp1(iREG,jREG) ) - zglamvGLO( zziv(iREG,jREG) , zzjvm1(iREG,jREG)) ) &
  &                * cos( gphivGLO( zziv(iREG,jREG) , zzjv(iREG,jREG) ) * rad )                                                )

  !-- local angle between the two grids :
  angleXv(iREG,jREG) = angoutXv(iREG,jREG) - anginX
  angleYv(iREG,jREG) = angoutYt(iREG,jREG) - anginY
  if ( abs(angleXv(iREG,jREG)/rad) .gt. 45.0 .or. abs(angleYv(iREG,jREG)/rad) .gt. 45.0 ) then
    write(*,*) '@@@@@@ WARNING : angle between 2 grids > 45 deg at point of the bdy :', iREG, jREG
  endif     

  !------------------------------------------------------------
  if ( zglamvREG(iREG,jREG) .ge. zglamvGLO(zziv(iREG,jREG),zzjt(iREG,jREG)) )  then
    ziWv(iREG,jREG)   = zziv  (iREG,jREG) ! closest point westward
    ziEv(iREG,jREG)   = zzivp1(iREG,jREG) ! closest point eastward
  else
    ziWv(iREG,jREG)   = zzivm1(iREG,jREG) ! closest point westward
    ziEv(iREG,jREG)   = zziv  (iREG,jREG) ! closest point eastward
  endif
  !--
  if ( gphivREG(iREG,jREG) .ge. gphivGLO(zziv(iREG,jREG),zzjt(iREG,jREG)) )  then
    zjSv(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point southward
    zjNv(iREG,jREG)   = zzjtp1(iREG,jREG) ! closest point northward
  else
    zjSv(iREG,jREG)   = zzjtm1(iREG,jREG) ! closest point southward
    zjNv(iREG,jREG)   = zzjt  (iREG,jREG) ! closest point northward
  endif
  !--
  distSW = sqrt(   ( cos(rad*gphivREG(iREG,jREG)) * ( zglamvREG(iREG,jREG) - zglamvGLO(ziWv(iREG,jREG),zjSv(iREG,jREG)) ) )**2  &
      &          + (                                     gphivREG(iREG,jREG) -  gphivGLO(ziWv(iREG,jREG),zjSv(iREG,jREG)) )**2  )
  distSE = sqrt(   ( cos(rad*gphivREG(iREG,jREG)) * ( zglamvREG(iREG,jREG) - zglamvGLO(ziEv(iREG,jREG),zjSv(iREG,jREG)) ) )**2  &
      &          + (                                     gphivREG(iREG,jREG) -  gphivGLO(ziEv(iREG,jREG),zjSv(iREG,jREG)) )**2  )
  distNW = sqrt(   ( cos(rad*gphivREG(iREG,jREG)) * ( zglamvREG(iREG,jREG) - zglamvGLO(ziWv(iREG,jREG),zjNv(iREG,jREG)) ) )**2  &
      &          + (                                     gphivREG(iREG,jREG) -  gphivGLO(ziWv(iREG,jREG),zjNv(iREG,jREG)) )**2  )
  distNE = sqrt(   ( cos(rad*gphivREG(iREG,jREG)) * ( zglamvREG(iREG,jREG) - zglamvGLO(ziEv(iREG,jREG),zjNv(iREG,jREG)) ) )**2  &
      &          + (                                     gphivREG(iREG,jREG) -  gphivGLO(ziEv(iREG,jREG),zjNv(iREG,jREG)) )**2  )
  !--
  wgSWv(iREG,jREG) =   ( zglamvGLO(ziEv(iREG,jREG),zjNv(iREG,jREG)) - zglamvREG(iREG,jREG) ) &
  &                  * (  gphivGLO(ziEv(iREG,jREG),zjNv(iREG,jREG)) -  gphivREG(iREG,jREG) )
  wgSEv(iREG,jREG) =   ( zglamvGLO(ziWv(iREG,jREG),zjNv(iREG,jREG)) - zglamvREG(iREG,jREG) ) &
  &                  * (  gphivREG(iREG,jREG) -  gphivGLO(ziWv(iREG,jREG),zjNv(iREG,jREG)) )
  wgNWv(iREG,jREG) =   ( zglamvREG(iREG,jREG) - zglamvGLO(ziEv(iREG,jREG),zjSv(iREG,jREG)) ) &
  &                  * (  gphivGLO(ziEv(iREG,jREG),zjSv(iREG,jREG)) -  gphivREG(iREG,jREG) )
  wgNEv(iREG,jREG) =   ( zglamvREG(iREG,jREG) - zglamvGLO(ziWv(iREG,jREG),zjSv(iREG,jREG)) ) &
  &                  * (  gphivREG(iREG,jREG) -  gphivGLO(ziWv(iREG,jREG),zjSv(iREG,jREG)) )

ENDDO !-- jREG
ENDDO !-- iREG

!==================================================================================
! 4- WRITE EVERYTHING NEEDED FOR INTERPOLATION IN A NETCDF FILE
!==================================================================================

write(*,*) 'Writing ', TRIM(file_coeff)

!status = NF90_CREATE(TRIM(file_coeff),or(NF90_NOCLOBBER,NF90_64BIT_OFFSET),fidcoeff)
status = NF90_CREATE(TRIM(file_coeff),NF90_NOCLOBBER,fidcoeff)
call erreur(status,.TRUE.,'create_weights_file')

!-- File dimensions 
status = NF90_DEF_DIM(fidcoeff,"x",mx_REG,dimID_x) ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidcoeff,"y",my_REG,dimID_y) ; call erreur(status,.TRUE.,"def_dimID_y")
                  
!-- Variables definition                            
status = NF90_DEF_VAR(fidcoeff,"glamt",NF90_FLOAT,(/dimID_x,dimID_y/),glamtREG_ID)  ; call erreur(status,.TRUE.,"def_var_glamtREG_ID")
status = NF90_DEF_VAR(fidcoeff,"gphit",NF90_FLOAT,(/dimID_x,dimID_y/),gphitREG_ID)  ; call erreur(status,.TRUE.,"def_var_gphitREG_ID")
status = NF90_DEF_VAR(fidcoeff,"glamf",NF90_FLOAT,(/dimID_x,dimID_y/),glamfREG_ID)  ; call erreur(status,.TRUE.,"def_var_glamfREG_ID")
status = NF90_DEF_VAR(fidcoeff,"gphif",NF90_FLOAT,(/dimID_x,dimID_y/),gphifREG_ID)  ; call erreur(status,.TRUE.,"def_var_gphifREG_ID")
status = NF90_DEF_VAR(fidcoeff,"glamu",NF90_FLOAT,(/dimID_x,dimID_y/),glamuREG_ID)  ; call erreur(status,.TRUE.,"def_var_glamuREG_ID")
status = NF90_DEF_VAR(fidcoeff,"gphiu",NF90_FLOAT,(/dimID_x,dimID_y/),gphiuREG_ID)  ; call erreur(status,.TRUE.,"def_var_gphiuREG_ID")
status = NF90_DEF_VAR(fidcoeff,"glamv",NF90_FLOAT,(/dimID_x,dimID_y/),glamvREG_ID)  ; call erreur(status,.TRUE.,"def_var_glamvREG_ID")
status = NF90_DEF_VAR(fidcoeff,"gphiv",NF90_FLOAT,(/dimID_x,dimID_y/),gphivREG_ID)  ; call erreur(status,.TRUE.,"def_var_gphivREG_ID")
!-
status = NF90_DEF_VAR(fidcoeff,"ziWt", NF90_INT  ,(/dimID_x,dimID_y/),ziWt_ID)    ; call erreur(status,.TRUE.,"def_var_ziWt_ID")
status = NF90_DEF_VAR(fidcoeff,"ziEt", NF90_INT  ,(/dimID_x,dimID_y/),ziEt_ID)    ; call erreur(status,.TRUE.,"def_var_ziEt_ID")
status = NF90_DEF_VAR(fidcoeff,"zjSt", NF90_INT  ,(/dimID_x,dimID_y/),zjSt_ID)    ; call erreur(status,.TRUE.,"def_var_zjSt_ID")
status = NF90_DEF_VAR(fidcoeff,"zjNt", NF90_INT  ,(/dimID_x,dimID_y/),zjNt_ID)    ; call erreur(status,.TRUE.,"def_var_zjNt_ID")
status = NF90_DEF_VAR(fidcoeff,"angXt",NF90_FLOAT,(/dimID_x,dimID_y/),angleXt_ID) ; call erreur(status,.TRUE.,"def_var_angleXt_ID")
status = NF90_DEF_VAR(fidcoeff,"angYt",NF90_FLOAT,(/dimID_x,dimID_y/),angleYt_ID) ; call erreur(status,.TRUE.,"def_var_angleYt_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNEt", NF90_FLOAT,(/dimID_x,dimID_y/),wgNEt_ID) ; call erreur(status,.TRUE.,"def_var_wgNEt_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSEt", NF90_FLOAT,(/dimID_x,dimID_y/),wgSEt_ID) ; call erreur(status,.TRUE.,"def_var_wgSEt_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSWt", NF90_FLOAT,(/dimID_x,dimID_y/),wgSWt_ID) ; call erreur(status,.TRUE.,"def_var_wgSWt_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNWt", NF90_FLOAT,(/dimID_x,dimID_y/),wgNWt_ID) ; call erreur(status,.TRUE.,"def_var_wgNWt_ID")
!-
status = NF90_DEF_VAR(fidcoeff,"ziWf", NF90_INT  ,(/dimID_x,dimID_y/),ziWf_ID)    ; call erreur(status,.TRUE.,"def_var_ziWf_ID")
status = NF90_DEF_VAR(fidcoeff,"ziEf", NF90_INT  ,(/dimID_x,dimID_y/),ziEf_ID)    ; call erreur(status,.TRUE.,"def_var_ziEf_ID")
status = NF90_DEF_VAR(fidcoeff,"zjSf", NF90_INT  ,(/dimID_x,dimID_y/),zjSf_ID)    ; call erreur(status,.TRUE.,"def_var_zjSf_ID")
status = NF90_DEF_VAR(fidcoeff,"zjNf", NF90_INT  ,(/dimID_x,dimID_y/),zjNf_ID)    ; call erreur(status,.TRUE.,"def_var_zjNf_ID")
status = NF90_DEF_VAR(fidcoeff,"angXf",NF90_FLOAT,(/dimID_x,dimID_y/),angleXf_ID) ; call erreur(status,.TRUE.,"def_var_angleXf_ID")
status = NF90_DEF_VAR(fidcoeff,"angYf",NF90_FLOAT,(/dimID_x,dimID_y/),angleYf_ID) ; call erreur(status,.TRUE.,"def_var_angleYf_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNEf", NF90_FLOAT,(/dimID_x,dimID_y/),wgNEf_ID) ; call erreur(status,.TRUE.,"def_var_wgNEf_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSEf", NF90_FLOAT,(/dimID_x,dimID_y/),wgSEf_ID) ; call erreur(status,.TRUE.,"def_var_wgSEf_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSWf", NF90_FLOAT,(/dimID_x,dimID_y/),wgSWf_ID) ; call erreur(status,.TRUE.,"def_var_wgSWf_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNWf", NF90_FLOAT,(/dimID_x,dimID_y/),wgNWf_ID) ; call erreur(status,.TRUE.,"def_var_wgNWf_ID")
!-
status = NF90_DEF_VAR(fidcoeff,"ziWu", NF90_INT  ,(/dimID_x,dimID_y/),ziWu_ID)    ; call erreur(status,.TRUE.,"def_var_ziWu_ID")
status = NF90_DEF_VAR(fidcoeff,"ziEu", NF90_INT  ,(/dimID_x,dimID_y/),ziEu_ID)    ; call erreur(status,.TRUE.,"def_var_ziEu_ID")
status = NF90_DEF_VAR(fidcoeff,"zjSu", NF90_INT  ,(/dimID_x,dimID_y/),zjSu_ID)    ; call erreur(status,.TRUE.,"def_var_zjSu_ID")
status = NF90_DEF_VAR(fidcoeff,"zjNu", NF90_INT  ,(/dimID_x,dimID_y/),zjNu_ID)    ; call erreur(status,.TRUE.,"def_var_zjNu_ID")
status = NF90_DEF_VAR(fidcoeff,"angXu",NF90_FLOAT,(/dimID_x,dimID_y/),angleXu_ID) ; call erreur(status,.TRUE.,"def_var_angleXu_ID")
status = NF90_DEF_VAR(fidcoeff,"angYu",NF90_FLOAT,(/dimID_x,dimID_y/),angleYu_ID) ; call erreur(status,.TRUE.,"def_var_angleYu_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNEu", NF90_FLOAT,(/dimID_x,dimID_y/),wgNEu_ID) ; call erreur(status,.TRUE.,"def_var_wgNEu_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSEu", NF90_FLOAT,(/dimID_x,dimID_y/),wgSEu_ID) ; call erreur(status,.TRUE.,"def_var_wgSEu_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSWu", NF90_FLOAT,(/dimID_x,dimID_y/),wgSWu_ID) ; call erreur(status,.TRUE.,"def_var_wgSWu_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNWu", NF90_FLOAT,(/dimID_x,dimID_y/),wgNWu_ID) ; call erreur(status,.TRUE.,"def_var_wgNWu_ID")
!-
status = NF90_DEF_VAR(fidcoeff,"ziWv", NF90_INT  ,(/dimID_x,dimID_y/),ziWv_ID)    ; call erreur(status,.TRUE.,"def_var_ziWv_ID")
status = NF90_DEF_VAR(fidcoeff,"ziEv", NF90_INT  ,(/dimID_x,dimID_y/),ziEv_ID)    ; call erreur(status,.TRUE.,"def_var_ziEv_ID")
status = NF90_DEF_VAR(fidcoeff,"zjSv", NF90_INT  ,(/dimID_x,dimID_y/),zjSv_ID)    ; call erreur(status,.TRUE.,"def_var_zjSv_ID")
status = NF90_DEF_VAR(fidcoeff,"zjNv", NF90_INT  ,(/dimID_x,dimID_y/),zjNv_ID)    ; call erreur(status,.TRUE.,"def_var_zjNv_ID")
status = NF90_DEF_VAR(fidcoeff,"angXv",NF90_FLOAT,(/dimID_x,dimID_y/),angleXv_ID) ; call erreur(status,.TRUE.,"def_var_angleXv_ID")
status = NF90_DEF_VAR(fidcoeff,"angYv",NF90_FLOAT,(/dimID_x,dimID_y/),angleYv_ID) ; call erreur(status,.TRUE.,"def_var_angleYv_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNEv", NF90_FLOAT,(/dimID_x,dimID_y/),wgNEv_ID) ; call erreur(status,.TRUE.,"def_var_wgNEv_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSEv", NF90_FLOAT,(/dimID_x,dimID_y/),wgSEv_ID) ; call erreur(status,.TRUE.,"def_var_wgSEv_ID")
status = NF90_DEF_VAR(fidcoeff,"wgSWv", NF90_FLOAT,(/dimID_x,dimID_y/),wgSWv_ID) ; call erreur(status,.TRUE.,"def_var_wgSWv_ID")
status = NF90_DEF_VAR(fidcoeff,"wgNWv", NF90_FLOAT,(/dimID_x,dimID_y/),wgNWv_ID) ; call erreur(status,.TRUE.,"def_var_wgNWv_ID")

!-- Global attribute
status = NF90_PUT_ATT(fidcoeff,NF90_GLOBAL,"history","built using build_linear_interp_weights.f90")
status = NF90_PUT_ATT(fidcoeff,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2")
call erreur(status,.TRUE.,"att_global")

!-- End of definitions
status = NF90_ENDDEF(fidcoeff) ; call erreur(status,.TRUE.,"end_definition") 


!-- Values to put in each variable 
status = NF90_PUT_VAR(fidcoeff,glamtREG_ID,glamtREG) ; call erreur(status,.TRUE.,"var_glamtREG_ID")
status = NF90_PUT_VAR(fidcoeff,gphitREG_ID,gphitREG) ; call erreur(status,.TRUE.,"var_gphitREG_ID")
status = NF90_PUT_VAR(fidcoeff,glamfREG_ID,glamfREG) ; call erreur(status,.TRUE.,"var_glamfREG_ID")
status = NF90_PUT_VAR(fidcoeff,gphifREG_ID,gphifREG) ; call erreur(status,.TRUE.,"var_gphifREG_ID")
status = NF90_PUT_VAR(fidcoeff,glamuREG_ID,glamuREG) ; call erreur(status,.TRUE.,"var_glamuREG_ID")
status = NF90_PUT_VAR(fidcoeff,gphiuREG_ID,gphiuREG) ; call erreur(status,.TRUE.,"var_gphiuREG_ID")
status = NF90_PUT_VAR(fidcoeff,glamvREG_ID,glamvREG) ; call erreur(status,.TRUE.,"var_glamvREG_ID")
status = NF90_PUT_VAR(fidcoeff,gphivREG_ID,gphivREG) ; call erreur(status,.TRUE.,"var_gphivREG_ID")
!-
status = NF90_PUT_VAR(fidcoeff,ziWt_ID,ziWt)         ; call erreur(status,.TRUE.,"var_ziWt_ID")
status = NF90_PUT_VAR(fidcoeff,ziEt_ID,ziEt)         ; call erreur(status,.TRUE.,"var_ziEt_ID")
status = NF90_PUT_VAR(fidcoeff,zjSt_ID,zjSt)         ; call erreur(status,.TRUE.,"var_zjSt_ID")
status = NF90_PUT_VAR(fidcoeff,zjNt_ID,zjNt)         ; call erreur(status,.TRUE.,"var_zjNt_ID")
status = NF90_PUT_VAR(fidcoeff,wgNEt_ID,wgNEt)       ; call erreur(status,.TRUE.,"var_wgNEt_ID")
status = NF90_PUT_VAR(fidcoeff,wgSEt_ID,wgSEt)       ; call erreur(status,.TRUE.,"var_wgSEt_ID")
status = NF90_PUT_VAR(fidcoeff,wgSWt_ID,wgSWt)       ; call erreur(status,.TRUE.,"var_wgSWt_ID")
status = NF90_PUT_VAR(fidcoeff,wgNWt_ID,wgNWt)       ; call erreur(status,.TRUE.,"var_wgNWt_ID")
status = NF90_PUT_VAR(fidcoeff,angleXt_ID,angleXt)   ; call erreur(status,.TRUE.,"var_angleXt_ID")
status = NF90_PUT_VAR(fidcoeff,angleYt_ID,angleYt)   ; call erreur(status,.TRUE.,"var_angleYt_ID")
!-
status = NF90_PUT_VAR(fidcoeff,ziWf_ID,ziWf)         ; call erreur(status,.TRUE.,"var_ziWf_ID")
status = NF90_PUT_VAR(fidcoeff,ziEf_ID,ziEf)         ; call erreur(status,.TRUE.,"var_ziEf_ID")
status = NF90_PUT_VAR(fidcoeff,zjSf_ID,zjSf)         ; call erreur(status,.TRUE.,"var_zjSf_ID")
status = NF90_PUT_VAR(fidcoeff,zjNf_ID,zjNf)         ; call erreur(status,.TRUE.,"var_zjNf_ID")
status = NF90_PUT_VAR(fidcoeff,wgNEf_ID,wgNEf)       ; call erreur(status,.TRUE.,"var_wgNEf_ID")
status = NF90_PUT_VAR(fidcoeff,wgSEf_ID,wgSEf)       ; call erreur(status,.TRUE.,"var_wgSEf_ID")
status = NF90_PUT_VAR(fidcoeff,wgSWf_ID,wgSWf)       ; call erreur(status,.TRUE.,"var_wgSWf_ID")
status = NF90_PUT_VAR(fidcoeff,wgNWf_ID,wgNWf)       ; call erreur(status,.TRUE.,"var_wgNWf_ID")
status = NF90_PUT_VAR(fidcoeff,angleXf_ID,angleXf)   ; call erreur(status,.TRUE.,"var_angleXf_ID")
status = NF90_PUT_VAR(fidcoeff,angleYf_ID,angleYf)   ; call erreur(status,.TRUE.,"var_angleYf_ID")
!-
status = NF90_PUT_VAR(fidcoeff,ziWu_ID,ziWu)         ; call erreur(status,.TRUE.,"var_ziWu_ID")
status = NF90_PUT_VAR(fidcoeff,ziEu_ID,ziEu)         ; call erreur(status,.TRUE.,"var_ziEu_ID")
status = NF90_PUT_VAR(fidcoeff,zjSu_ID,zjSu)         ; call erreur(status,.TRUE.,"var_zjSu_ID")
status = NF90_PUT_VAR(fidcoeff,zjNu_ID,zjNu)         ; call erreur(status,.TRUE.,"var_zjNu_ID")
status = NF90_PUT_VAR(fidcoeff,wgNEu_ID,wgNEu)       ; call erreur(status,.TRUE.,"var_wgNEu_ID")
status = NF90_PUT_VAR(fidcoeff,wgSEu_ID,wgSEu)       ; call erreur(status,.TRUE.,"var_wgSEu_ID")
status = NF90_PUT_VAR(fidcoeff,wgSWu_ID,wgSWu)       ; call erreur(status,.TRUE.,"var_wgSWu_ID")
status = NF90_PUT_VAR(fidcoeff,wgNWu_ID,wgNWu)       ; call erreur(status,.TRUE.,"var_wgNWu_ID")
status = NF90_PUT_VAR(fidcoeff,angleXu_ID,angleXu)   ; call erreur(status,.TRUE.,"var_angleXu_ID")
status = NF90_PUT_VAR(fidcoeff,angleYu_ID,angleYu)   ; call erreur(status,.TRUE.,"var_angleYu_ID")
!-
status = NF90_PUT_VAR(fidcoeff,ziWv_ID,ziWv)         ; call erreur(status,.TRUE.,"var_ziWv_ID")
status = NF90_PUT_VAR(fidcoeff,ziEv_ID,ziEv)         ; call erreur(status,.TRUE.,"var_ziEv_ID")
status = NF90_PUT_VAR(fidcoeff,zjSv_ID,zjSv)         ; call erreur(status,.TRUE.,"var_zjSv_ID")
status = NF90_PUT_VAR(fidcoeff,zjNv_ID,zjNv)         ; call erreur(status,.TRUE.,"var_zjNv_ID")
status = NF90_PUT_VAR(fidcoeff,wgNEv_ID,wgNEv)       ; call erreur(status,.TRUE.,"var_wgNEv_ID")
status = NF90_PUT_VAR(fidcoeff,wgSEv_ID,wgSEv)       ; call erreur(status,.TRUE.,"var_wgSEv_ID")
status = NF90_PUT_VAR(fidcoeff,wgSWv_ID,wgSWv)       ; call erreur(status,.TRUE.,"var_wgSWv_ID")
status = NF90_PUT_VAR(fidcoeff,wgNWv_ID,wgNWv)       ; call erreur(status,.TRUE.,"var_wgNWv_ID")
status = NF90_PUT_VAR(fidcoeff,angleXv_ID,angleXv)   ; call erreur(status,.TRUE.,"var_angleXv_ID")
status = NF90_PUT_VAR(fidcoeff,angleYv_ID,angleYv)   ; call erreur(status,.TRUE.,"var_angleYv_ID")

!-- End of writing                           
status = NF90_CLOSE(fidcoeff) ; call erreur(status,.TRUE.,"end_writing_weights")

end program buildweight


!==================================================================================


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
