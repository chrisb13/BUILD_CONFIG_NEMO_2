program modif                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, LGGE-CNRS, Feb. 2015
!
! purpose: extract iceberg runoff from ORCA025 grid to AMU12
!
! 0- Initializations
! 1- Read information on grids
! 2- Read iceberg runoff
! 3- Read REGIONAL mesh/mask file
! 4- Read mesh_mask of global grid
! 5- Projection onto regional grid
! 6- Define socoefr to avoid SSS restoring within some distance from the coast
! 7- Writing initial state for salinity
!
! history: - March 2017: GitHub version with namelist (N. Jourdain, IGE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /runoff/ file_runoff_in, file_mask_runoff, nn_band
INTEGER                                  :: nn_band
CHARACTER(LEN=50)                        :: config
CHARACTER(LEN=150)                       :: file_runoff_in, file_mask_runoff, config_dir

INTEGER                                  :: fidA, status, dimID_x, dimID_y, dimID_time_counter, dimID_z, dimID_t, mx, my, mtime_counter, &
&                                           mzreg, mtreg, time_counter_ID, nav_lon_ID, nav_lat_ID, Melt_ID, fidM, fidglo, jmin_ORCA12,   &
&                                           i, j, iGLO, jGLO, kk, rr, rs, mxglo, myglo, mzglo, mtglo, fidMSH, l, tmask_ID, mb, kki, kkj, &
&                                           socoefr_ID, fidC, ai, aj, bi, bj, iREG, jREG, mx_REG, my_REG, imin_ORCA12, dij, im1, ip1,    &
&                                           jm1, jp1, mx_tmp, my_tmp
CHARACTER(LEN=150)                       :: file_runoff_out, file_in_mask_REG, file_in_coord_REG
REAL*8                                   :: chkland
REAL*4,ALLOCATABLE,DIMENSION(:)          :: time_counter           
REAL*4,ALLOCATABLE,DIMENSION(:,:)        :: nav_lon, nav_lat, lonreg, latreg, socoefr
REAL*4,ALLOCATABLE,DIMENSION(:,:,:)      :: Melt, Meltreg, tmp_Meltreg
INTEGER*1,ALLOCATABLE,DIMENSION(:,:)     :: tmask_REG, tmask_GLO
INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmask

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

! Default values (replaced with namelist values if specified):
config_dir        = '.'

!- read namelist values :
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=runoff)
CLOSE(1)

!- name of regional mesh_mask (input) :
write(file_in_mask_REG,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/mesh_mask_',a,'.nc')

!- name of regional coordinates (input) :
write(file_in_coord_REG,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

! name of interpolation weights file :
write(file_coeff,103) TRIM(config_dir), TRIM(config)
103 FORMAT(a,'/coeff_linear_',a,'.nc')

!- output file names :
write(file_runoff_out,201)  TRIM(config_dir), TRIM(config)
201 FORMAT(a,'/runoff_iceberg_',a,'.nc')

!=================================================================================
! 1- Read iceberg runoff :
!=================================================================================

write(*,*) 'Reading ', TRIM(file_runoff_in)
 
status = NF90_OPEN(TRIM(file_runoff_in),0,fidA)  ; call erreur(status,.TRUE.,"read input runoff") 
                                   
status = NF90_INQ_DIMID(fidA,"x",dimID_x)                        ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"y",dimID_y)                        ; call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"time_counter",dimID_time_counter)  ; call erreur(status,.TRUE.,"inq_dimID_time_counter")
                                       
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx)                       ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my)                       ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_time_counter,len=mtime_counter) ; call erreur(status,.TRUE.,"inq_dim_time_counter")
       
ALLOCATE(  time_counter(mtime_counter)  ) 
ALLOCATE(  nav_lon(mx,my)  ) 
ALLOCATE(  nav_lat(mx,my)  ) 
ALLOCATE(  Melt(mx,my,mtime_counter)  ) 
         
status = NF90_INQ_VARID(fidA,"time_counter",time_counter_ID) ; call erreur(status,.TRUE.,"inq_time_counter_ID")
status = NF90_INQ_VARID(fidA,"nav_lon",nav_lon_ID)           ; call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)           ; call erreur(status,.TRUE.,"inq_nav_lat_ID")
status = NF90_INQ_VARID(fidA,"Melt",Melt_ID)                 ; call erreur(status,.TRUE.,"inq_Melt_ID")
                                      
status = NF90_GET_VAR(fidA,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"getvar_time_counter")
status = NF90_GET_VAR(fidA,nav_lon_ID,nav_lon)           ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidA,nav_lat_ID,nav_lat)           ; call erreur(status,.TRUE.,"getvar_nav_lat")
status = NF90_GET_VAR(fidA,Melt_ID,Melt)                 ; call erreur(status,.TRUE.,"getvar_Melt")
                              
status = NF90_CLOSE(fidA) ; call erreur(status,.TRUE.,"fin_lecture")     
                                      
!=================================================================================
! 2- Read REGIONAL mesh/mask file
!=================================================================================

write(*,*) 'Reading ', TRIM(file_in_mask_REG)

status = NF90_OPEN(TRIM(file_in_mask_REG),0,fidMSH) ; call erreur(status,.TRUE.,"read AMU12 mesh mask") 

status = NF90_INQ_DIMID(fidMSH,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidMSH,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidMSH,"z",dimID_z) ; call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidMSH,"t",dimID_t) ; call erreur(status,.TRUE.,"inq_dimID_t")

status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_x,len=mx_REG) ; call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_y,len=my_REG) ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_z,len=mzreg) ; call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_t,len=mtreg) ; call erreur(status,.TRUE.,"inq_dim_t")

ALLOCATE(  tmask(mx_REG,my_REG,mzreg,mtreg)  )
ALLOCATE(  tmask_REG(mx_REG,my_REG)  )
ALLOCATE(  Meltreg(mx_REG,my_REG,12) )
ALLOCATE(  tmp_Meltreg(mx_REG,my_REG,12) )
ALLOCATE(  lonreg(mx_REG,my_REG)  ) 
ALLOCATE(  latreg(mx_REG,my_REG)  ) 

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)     ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lon",nav_lon_ID) ; call erreur(status,.TRUE.,"inq_nav_lon_ID")
status = NF90_INQ_VARID(fidMSH,"nav_lat",nav_lat_ID) ; call erreur(status,.TRUE.,"inq_nav_lat_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)         ; call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,nav_lon_ID,lonreg)      ; call erreur(status,.TRUE.,"getvar_nav_lon")
status = NF90_GET_VAR(fidMSH,nav_lat_ID,latreg)      ; call erreur(status,.TRUE.,"getvar_nav_lat")

tmask_REG(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

status = NF90_CLOSE(fidMSH) ; call erreur(status,.TRUE.,"fin_lecture")     

!=================================================================================
! 3- Read mesh_mask of global grid
!=================================================================================

write(*,*) 'Reading ', TRIM(file_mask_runoff)
                                           
status = NF90_OPEN(TRIM(file_mask_runoff),0,fidglo) ; call erreur(status,.TRUE.,"read large-scale mesh") 
                                           
status = NF90_INQ_DIMID(fidglo,"t",dimID_t) ; call erreur(status,.TRUE.,"inq_dimID_t")
status = NF90_INQ_DIMID(fidglo,"z",dimID_z) ; call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidglo,"y",dimID_y) ; call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidglo,"x",dimID_x) ; call erreur(status,.TRUE.,"inq_dimID_x")
                                               
status = NF90_INQUIRE_DIMENSION(fidglo,dimID_t,len=mtglo) ; call erreur(status,.TRUE.,"inq_dim_t")
status = NF90_INQUIRE_DIMENSION(fidglo,dimID_z,len=mzglo) ; call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidglo,dimID_y,len=myglo) ; call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidglo,dimID_x,len=mxglo) ; call erreur(status,.TRUE.,"inq_dim_x")
               
ALLOCATE(  tmask(mxglo,myglo,mzglo,mtglo)  ) 
ALLOCATE(  tmask_GLO(mxglo,myglo)  ) 
status = NF90_INQ_VARID(fidglo,"tmask",tmask_ID)  ; call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_GET_VAR(fidglo,tmask_ID,tmask)      ; call erreur(status,.TRUE.,"getvar_tmask")
                                      
status = NF90_CLOSE(fidglo) ; call erreur(status,.TRUE.,"fin_lecture")     

if ( mx .ne. mxglo .or. my.ne.myglo ) then
  write(*,*) 'PROBLEM WITH DIMENSIONS >>>>>>> stop !'
  stop
endif

tmask_GLO(:,:)=tmask(:,:,1,1)
DEALLOCATE(tmask)

!=================================================================================
! 4- Read interpolation weights
!=================================================================================

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
! 5- Projection onto regional grid :
!=================================================================================


do iREG=1,mx_REG
do jREG=1,my_REG
  aNW = tmask_GLO( ziWt(iREG,jREG), zjNt(iREG,jREG) ) * wgNWt(iREG,jREG)
  aSW = tmask_GLO( ziWt(iREG,jREG), zjSt(iREG,jREG) ) * wgSWt(iREG,jREG)
  aNE = tmask_GLO( ziEt(iREG,jREG), zjNt(iREG,jREG) ) * wgNEt(iREG,jREG)
  aSE = tmask_GLO( ziEt(iREG,jREG), zjSt(iREG,jREG) ) * wgSEt(iREG,jREG)
  aa = aNW+aSW+aNE+aSE
  if ( aa .gt. eps .and. zjSt(iREG,jREG) .gt. 1 ) then
    do l=1,mtime_counter
      Meltreg(iREG,jREG,l)  =  (   Melt( ziWt(iREG,jREG), zjNt(iREG,jREG), l ) * aNW   &
      &                          + Melt( ziWt(iREG,jREG), zjSt(iREG,jREG), l ) * aSW   &
      &                          + Melt( ziEt(iREG,jREG), zjNt(iREG,jREG), l ) * aNE   &
      &                          + Melt( ziEt(iREG,jREG), zjSt(iREG,jREG), l ) * aSE ) &
      &                        * tmask_REG(iREG,jREG) / aa
    enddo
  else
    Meltreg(iREG,jREG,:)  = 0.e0
  endif
enddo
enddo

!=================================================================================
! 6- Define socoefr to avoid SSS restoring within some distance from the coast:
!    socoefr=0.5 where no SSS relaxation is applied, and = 0.0 everywhere else.
!    (and between 0 and 0.5 for partial SSS relaxation).
!=================================================================================

mb=nn_band

ALLOCATE( socoefr(mx_REG,my_REG) )
socoefr(:,:) = 0.0

!! manual corrections :
if ( TRIM(config) == 'AMU12' .or. TRIM(config) == 'AMU12y' ) then
  tmask_REG(240:320,1:80) = 0.0 ! to remove SSS restoring in the TG & PIG Bay
elseif ( TRIM(config) == 'WED12' ) then
  tmask_REG(:,1:415) = 0.0 ! region where the old ORCA025 grid was masked
endif

DO i=1,mx_REG
DO j=1,my_REG
  ! mask surrounding domain excluded :
  chkland=SUM(SUM(1.00000000000000*(1-tmask_REG(MAX(2,i-mb):MIN(i+mb,mx_REG-1),MAX(2,j-mb):MIN(j+mb,my_REG-1))),2),1)
  if ( chkland .gt. 5.5 ) socoefr(i,j) = 0.5 * MIN(chkland,0.15*4*mb**2) / (0.15*4*mb**2)  ! 0.15 is for 15% of the square with land
ENDDO
ENDDO

DO i=1,mx_REG
DO j=1,my_REG
  if ( tmask_REG(i,j) .eq. 0 ) socoefr(i,j) = 0.5  ! Important to mask where there are ice shelf cavities !
ENDDO
ENDDO

!=================================================================================
! 7- Writing new runoff file on REGIONAL grid
!=================================================================================
                                      
status = NF90_CREATE(TRIM(file_runoff_out),NF90_NOCLOBBER,fidM) ; call erreur(status,.TRUE.,'create new runoff file')                     
                                        
status = NF90_DEF_DIM(fidM,"x",mx_REG,dimID_x)                               ; call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my_REG,dimID_y)                               ; call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"time_counter",NF90_UNLIMITED,dimID_time_counter) ; call erreur(status,.TRUE.,"def_dimID_time_counter")
                                      
status = NF90_DEF_VAR(fidM,"time_counter",NF90_FLOAT,(/dimID_time_counter/),time_counter_ID) ; call erreur(status,.TRUE.,"def_var_time_counter_ID")
status = NF90_DEF_VAR(fidM,"nav_lon",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lon_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lon_ID")
status = NF90_DEF_VAR(fidM,"nav_lat",NF90_FLOAT,(/dimID_x,dimID_y/),nav_lat_ID)              ; call erreur(status,.TRUE.,"def_var_nav_lat_ID")
status = NF90_DEF_VAR(fidM,"Melt",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),Melt_ID) ; call erreur(status,.TRUE.,"def_var_Melt_ID")
status = NF90_DEF_VAR(fidM,"socoefr",NF90_FLOAT,(/dimID_x,dimID_y/),socoefr_ID)              ; call erreur(status,.TRUE.,"def_var_socoefr_ID")
           
status = NF90_PUT_ATT(fidM,time_counter_ID,"long_name","Time axis")                     ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"title","Time")                              ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"time_origin","1989-01-01 00:00:00")         ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"units","seconds since 1989-01-01 00:00:00") ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,time_counter_ID,"calendar","noleap")                         ; call erreur(status,.TRUE.,"put_att_time_counter_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"nav_model","Default grid") ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"long_name","Longitude")    ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_max",180.)           ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"valid_min",-180.)          ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lon_ID,"units","degrees_east")     ; call erreur(status,.TRUE.,"put_att_nav_lon_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"nav_model","Default grid") ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"long_name","Latitude")     ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_max",90.)            ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"valid_min",-90.)           ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,nav_lat_ID,"units","degrees_north")    ; call erreur(status,.TRUE.,"put_att_nav_lat_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"savelog10",0.)                   ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"axis","T")                       ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"online_operation","N/A")         ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"short_name","Melt")              ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"long_name","Icebergs melt rate") ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,Melt_ID,"units","kg/m2/s")                ; call erreur(status,.TRUE.,"put_att_Melt_ID")
status = NF90_PUT_ATT(fidM,socoefr_ID,"units","-")            ; call erreur(status,.TRUE.,"put_att_socoefr_ID")
status = NF90_PUT_ATT(fidM,socoefr_ID,"short_name","socoefr") ; call erreur(status,.TRUE.,"put_att_socoefr_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using extract_runoff_icebergs.f90")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tools","https://github.com/nicojourdain/BUILD_CONFIG_NEMO")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
                              
status = NF90_ENDDEF(fidM) ; call erreur(status,.TRUE.,"fin_definition") 
                              
status = NF90_PUT_VAR(fidM,time_counter_ID,time_counter) ; call erreur(status,.TRUE.,"var_time_counter_ID")
status = NF90_PUT_VAR(fidM,nav_lon_ID,lonreg)            ; call erreur(status,.TRUE.,"var_nav_lon_ID")
status = NF90_PUT_VAR(fidM,nav_lat_ID,latreg)            ; call erreur(status,.TRUE.,"var_nav_lat_ID")
status = NF90_PUT_VAR(fidM,Melt_ID,Meltreg)              ; call erreur(status,.TRUE.,"var_Melt_ID")
status = NF90_PUT_VAR(fidM,socoefr_ID,socoefr)           ; call erreur(status,.TRUE.,"var_socoefr_ID")
                              
status = NF90_CLOSE(fidM) ;call erreur(status,.TRUE.,"final")         

end program modif

!----------------------------

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
