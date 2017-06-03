program modif                                         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, May 2017
!
! Script to extract the coordinates from a GLOBAL grid (e.g. eORCA12) to a
! REGIONAL grid (e.g. AMU12, WED12).
!
! 0- Initializations
! 1- Reading coordinates of global domain
! 2- Writing coordinates file for the regional domain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf                                            

IMPLICIT NONE                                         

!-- namelist parameters :
namelist /general/ config, config_dir
namelist /griddata/ inputdir, file_in_coord_extract, file_in_bathy_extract, ln_dateline, nn_perio,     &
&                   file_in_bathy_bdy, ln_isfcav, file_in_coord_bdy, nn_imin_extract, nn_imax_extract, &
&                   nn_jmin_extract, nn_jmax_extract
INTEGER                               :: nn_imin_extract, nn_imax_extract, nn_jmin_extract, nn_jmax_extract, nn_perio
CHARACTER(LEN=50)                     :: config
CHARACTER(LEN=150)                    :: inputdir, file_in_bathy_extract, file_in_coord_extract, file_in_bathy_bdy, &
&                                        config_dir, file_in_coord_bdy
LOGICAL                               :: ln_isfcav, ln_dateline

!-- local variables :
INTEGER :: fidORCA12, fidM, status, dimID_y, dimID_x, nav_lat_ID, nav_lon_ID, fidCOORDreg, & 
&          my_GLO, mx_GLO, my_REG, mx_REG, imin_GLO, imax_GLO, jmin_GLO, jmax_GLO,         &
&          iREG, jREG, jtmp, npts, kk, mx_tmp, my_tmp, e2f_ID, e2v_ID, fidGLO,             &
&          e2u_ID, e2t_ID, e1f_ID, e1v_ID, e1u_ID, e1t_ID, gphif_ID, gphiv_ID, gphiu_ID,   &
&          gphit_ID, glamf_ID, glamv_ID, glamu_ID, glamt_ID

CHARACTER(LEN=150) :: aaa, file_bathy_out, file_coord_out, command_str

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: nav_lat_ORCA12, nav_lon_ORCA12, nav_lat_ORCA025, nav_lon_ORCA025, nav_lat_REG, nav_lon_REG 

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: e2f_GLO, e2v_GLO, e2u_GLO, e2t_GLO, e1f_GLO, e1v_GLO, e1u_GLO, e1t_GLO, gphif_GLO,     &
&                                          gphiv_GLO, gphiu_GLO, gphit_GLO, glamf_GLO, glamv_GLO, glamu_GLO, glamt_GLO,           &
&                                          e2f_REG, e2v_REG, e2u_REG, e2t_REG, e1f_REG, e1v_REG, e1u_REG, e1t_REG,                &
&                                          gphif_REG, gphiv_REG, gphiu_REG, gphit_REG, glamf_REG, glamv_REG, glamu_REG, glamt_REG

INTEGER(KIND=4), DIMENSION(8) :: imin, imax, jmin

!=================================================================================
! 0- Initializations 
!=================================================================================

! Default values (replaced with namelist values if specified):
config_dir        = '.'

!- read namelist values
OPEN (UNIT=1, FILE='namelist_pre' )
READ (UNIT=1, NML=general)
READ (UNIT=1, NML=griddata)
CLOSE(1)

! name of regional bathymetry file (output file) :
write(file_bathy_out,101) TRIM(config_dir), TRIM(config)
101 FORMAT(a,'/bathy_meter_',a,'.nc')

! name of regional coordinates file (output file) :
write(file_coord_out,102) TRIM(config_dir), TRIM(config)
102 FORMAT(a,'/coordinates_',a,'.nc')

!-
imin_GLO = nn_imin_extract
imax_GLO = nn_imax_extract
jmin_GLO = nn_jmin_extract
jmax_GLO = nn_jmax_extract

mx_REG = imax_GLO - imin_GLO + 1
my_REG = jmax_GLO - jmin_GLO + 1

!- create the input directory is needed (config_dir)
write(command_str,888) TRIM(config_dir)
888 FORMAT('mkdir ',a)
CALL system(TRIM(command_str))

!=================================================================================
! 1- Reading coordinates on global domain
!=================================================================================

write(*,*) 'Extracting regional domain from : ', TRIM(file_in_coord_extract)

status = NF90_OPEN(TRIM(file_in_coord_extract),0,fidGLO); call erreur(status,.TRUE.,"read input coordinates") 

status = NF90_INQ_DIMID(fidGLO,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidGLO,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
  
status = NF90_INQUIRE_DIMENSION(fidGLO,dimID_y,len=my_GLO); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidGLO,dimID_x,len=mx_GLO); call erreur(status,.TRUE.,"inq_dim_x")

ALLOCATE(  e2f_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e2v_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e2u_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e2t_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e1f_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e1v_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e1u_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  e1t_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  gphif_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  gphiv_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  gphiu_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  gphit_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  glamf_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  glamv_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  glamu_GLO(mx_GLO,my_GLO)  ) 
ALLOCATE(  glamt_GLO(mx_GLO,my_GLO)  ) 

status = NF90_INQ_VARID(fidGLO,"e2f",e2f_ID);     call erreur(status,.TRUE.,"inq_e2f_ID")
status = NF90_INQ_VARID(fidGLO,"e2v",e2v_ID);     call erreur(status,.TRUE.,"inq_e2v_ID")
status = NF90_INQ_VARID(fidGLO,"e2u",e2u_ID);     call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidGLO,"e2t",e2t_ID);     call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidGLO,"e1f",e1f_ID);     call erreur(status,.TRUE.,"inq_e1f_ID")
status = NF90_INQ_VARID(fidGLO,"e1v",e1v_ID);     call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidGLO,"e1u",e1u_ID);     call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidGLO,"e1t",e1t_ID);     call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidGLO,"gphif",gphif_ID); call erreur(status,.TRUE.,"inq_gphif_ID")
status = NF90_INQ_VARID(fidGLO,"gphiv",gphiv_ID); call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidGLO,"gphiu",gphiu_ID); call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidGLO,"gphit",gphit_ID); call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidGLO,"glamf",glamf_ID); call erreur(status,.TRUE.,"inq_glamf_ID")
status = NF90_INQ_VARID(fidGLO,"glamv",glamv_ID); call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidGLO,"glamu",glamu_ID); call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidGLO,"glamt",glamt_ID); call erreur(status,.TRUE.,"inq_glamt_ID")
   
status = NF90_GET_VAR(fidGLO,e2f_ID,e2f_GLO);     call erreur(status,.TRUE.,"getvar_e2f")
status = NF90_GET_VAR(fidGLO,e2v_ID,e2v_GLO);     call erreur(status,.TRUE.,"getvar_e2v")
status = NF90_GET_VAR(fidGLO,e2u_ID,e2u_GLO);     call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidGLO,e2t_ID,e2t_GLO);     call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidGLO,e1f_ID,e1f_GLO);     call erreur(status,.TRUE.,"getvar_e1f")
status = NF90_GET_VAR(fidGLO,e1v_ID,e1v_GLO);     call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidGLO,e1u_ID,e1u_GLO);     call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidGLO,e1t_ID,e1t_GLO);     call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidGLO,gphif_ID,gphif_GLO); call erreur(status,.TRUE.,"getvar_gphif")
status = NF90_GET_VAR(fidGLO,gphiv_ID,gphiv_GLO); call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidGLO,gphiu_ID,gphiu_GLO); call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidGLO,gphit_ID,gphit_GLO); call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidGLO,glamf_ID,glamf_GLO); call erreur(status,.TRUE.,"getvar_glamf")
status = NF90_GET_VAR(fidGLO,glamv_ID,glamv_GLO); call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidGLO,glamu_ID,glamu_GLO); call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidGLO,glamt_ID,glamt_GLO); call erreur(status,.TRUE.,"getvar_glamt")

status = NF90_CLOSE(fidGLO); call erreur(status,.TRUE.,"End read coordinates")     
   
!=================================================================================
! 2- Writing coordinates file for the regional domain :                                   
!=================================================================================

ALLOCATE(  e2f_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2u_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e2t_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1f_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1v_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1u_REG(mx_REG,my_REG)  ) 
ALLOCATE(  e1t_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphif_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphiv_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphiu_REG(mx_REG,my_REG)  ) 
ALLOCATE(  gphit_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamf_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamv_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamu_REG(mx_REG,my_REG)  ) 
ALLOCATE(  glamt_REG(mx_REG,my_REG)  ) 

write(*,*) 'Creating ', TRIM(file_coord_out)
   
status = NF90_CREATE(TRIM(file_coord_out),NF90_NOCLOBBER,fidCOORDreg); call erreur(status,.TRUE.,'create file_coord_out') 

status = NF90_DEF_DIM(fidCOORDreg,"x",mx_REG,dimID_x); call erreur(status,.TRUE.,"def_dimID_coord_x")
status = NF90_DEF_DIM(fidCOORDreg,"y",my_REG,dimID_y); call erreur(status,.TRUE.,"def_dimID_coord_y")

status = NF90_DEF_VAR(fidCOORDreg,"e2f",NF90_DOUBLE,(/dimID_x,dimID_y/),e2f_ID);     call erreur(status,.TRUE.,"def_var_e2f_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2v",NF90_DOUBLE,(/dimID_x,dimID_y/),e2v_ID);     call erreur(status,.TRUE.,"def_var_e2v_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2u",NF90_DOUBLE,(/dimID_x,dimID_y/),e2u_ID);     call erreur(status,.TRUE.,"def_var_e2u_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e2t",NF90_DOUBLE,(/dimID_x,dimID_y/),e2t_ID);     call erreur(status,.TRUE.,"def_var_e2t_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1f",NF90_DOUBLE,(/dimID_x,dimID_y/),e1f_ID);     call erreur(status,.TRUE.,"def_var_e1f_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1v",NF90_DOUBLE,(/dimID_x,dimID_y/),e1v_ID);     call erreur(status,.TRUE.,"def_var_e1v_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1u",NF90_DOUBLE,(/dimID_x,dimID_y/),e1u_ID);     call erreur(status,.TRUE.,"def_var_e1u_ID")
status = NF90_DEF_VAR(fidCOORDreg,"e1t",NF90_DOUBLE,(/dimID_x,dimID_y/),e1t_ID);     call erreur(status,.TRUE.,"def_var_e1t_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphif",NF90_DOUBLE,(/dimID_x,dimID_y/),gphif_ID); call erreur(status,.TRUE.,"def_var_gphif_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphiv",NF90_DOUBLE,(/dimID_x,dimID_y/),gphiv_ID); call erreur(status,.TRUE.,"def_var_gphiv_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphiu",NF90_DOUBLE,(/dimID_x,dimID_y/),gphiu_ID); call erreur(status,.TRUE.,"def_var_gphiu_ID")
status = NF90_DEF_VAR(fidCOORDreg,"gphit",NF90_DOUBLE,(/dimID_x,dimID_y/),gphit_ID); call erreur(status,.TRUE.,"def_var_gphit_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamf",NF90_DOUBLE,(/dimID_x,dimID_y/),glamf_ID); call erreur(status,.TRUE.,"def_var_glamf_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamv",NF90_DOUBLE,(/dimID_x,dimID_y/),glamv_ID); call erreur(status,.TRUE.,"def_var_glamv_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamu",NF90_DOUBLE,(/dimID_x,dimID_y/),glamu_ID); call erreur(status,.TRUE.,"def_var_glamu_ID")
status = NF90_DEF_VAR(fidCOORDreg,"glamt",NF90_DOUBLE,(/dimID_x,dimID_y/),glamt_ID); call erreur(status,.TRUE.,"def_var_glamt_ID")

status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"history","Created using extract_bathy_coord.f90"); call erreur(status,.TRUE.,"put_att_GLOB1")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"domain",TRIM(aaa));                                call erreur(status,.TRUE.,"put_att_GLOB2")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imin_extraction",imin_GLO);                     call erreur(status,.TRUE.,"put_att_GLOB3")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"imax_extraction",imax_GLO);                     call erreur(status,.TRUE.,"put_att_GLOB4")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmin_extraction",jmin_GLO);                     call erreur(status,.TRUE.,"put_att_GLOB5")
status = NF90_PUT_ATT(fidCOORDreg,NF90_GLOBAL,"jmax_extraction",jmax_GLO);                     call erreur(status,.TRUE.,"put_att_GLOB6")

status = NF90_ENDDEF(fidCOORDreg); call erreur(status,.TRUE.,"end_definition_coord") 

e2f_REG   (:,:) = e2f_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e2v_REG   (:,:) = e2v_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e2u_REG   (:,:) = e2u_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e2t_REG   (:,:) = e2t_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e1f_REG   (:,:) = e1f_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e1v_REG   (:,:) = e1v_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e1u_REG   (:,:) = e1u_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
e1t_REG   (:,:) = e1t_GLO  (imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
gphif_REG (:,:) = gphif_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
gphiv_REG (:,:) = gphiv_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
gphiu_REG (:,:) = gphiu_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
gphit_REG (:,:) = gphit_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
glamf_REG (:,:) = glamf_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
glamv_REG (:,:) = glamv_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
glamu_REG (:,:) = glamu_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)
glamt_REG (:,:) = glamt_GLO(imin_GLO:imax_GLO,jmin_GLO:jmax_GLO)

status = NF90_PUT_VAR(fidCOORDreg,e2f_ID,e2f_REG);      call erreur(status,.TRUE.,"var_e2f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2v_ID,e2v_REG);      call erreur(status,.TRUE.,"var_e2v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2u_ID,e2u_REG);      call erreur(status,.TRUE.,"var_e2u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e2t_ID,e2t_REG);      call erreur(status,.TRUE.,"var_e2t_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1f_ID,e1f_REG);      call erreur(status,.TRUE.,"var_e1f_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1v_ID,e1v_REG);      call erreur(status,.TRUE.,"var_e1v_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1u_ID,e1u_REG);      call erreur(status,.TRUE.,"var_e1u_ID")
status = NF90_PUT_VAR(fidCOORDreg,e1t_ID,e1t_REG);      call erreur(status,.TRUE.,"var_e1t_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphif_ID,gphif_REG);  call erreur(status,.TRUE.,"var_gphif_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiv_ID,gphiv_REG);  call erreur(status,.TRUE.,"var_gphiv_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphiu_ID,gphiu_REG);  call erreur(status,.TRUE.,"var_gphiu_ID")
status = NF90_PUT_VAR(fidCOORDreg,gphit_ID,gphit_REG);  call erreur(status,.TRUE.,"var_gphit_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamf_ID,glamf_REG);  call erreur(status,.TRUE.,"var_glamf_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamv_ID,glamv_REG);  call erreur(status,.TRUE.,"var_glamv_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamu_ID,glamu_REG);  call erreur(status,.TRUE.,"var_glamu_ID")
status = NF90_PUT_VAR(fidCOORDreg,glamt_ID,glamt_REG);  call erreur(status,.TRUE.,"var_glamt_ID")

status = NF90_CLOSE(fidCOORDreg); call erreur(status,.TRUE.,"final")         

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
