#!/bin/bash
#####################################################################
#
# N. Jourdain, IGE-CNRS, Sep. 2017
#
# used to calculate the barotropic stream function (BSF).
#
# Before running this script, you need to compile the CDFTOOLS
# (see https://github.com/meom-group/CDFTOOLS)
# You may need to adapt the variable names in modcdfnames.F90
#
# IMPORTANT: choose OPTION to match with the considered simulation. 
#
#####################################################################

CDFTOOL_DIR='/home/njourd/CDFTOOLSv4'
CONFCASE="SOUTHORCA025.L75-GNM029"
INPUT_DIR="$SHAREDELMER/SOUTHORCA025_GNM029_TMP"
YEARi=1979
YEARf=1981
OPTION='-ssh'  # Either '-vvl'  (z*-coordinates, i.e. calculation based on time-varying e3[uv]) 
               #     or '-ssh'  (z-coordinate, with non-linear free surface, i.e. constant e3[uv] and ssh)
               #     or '    '  (z-coordinate, with linear free surf and constant e3[uv], not accounting for ssh)
               #     or '-full' (z-coordinate, full steps, with constant e3[uv], not accounting for ssh)

ln -s -v /store/njourd/SOUTHORCA025_GNM029/SOUTHORCA025_mesh_zgr.nc   mesh_zgr.nc
ln -s -v /store/njourd/SOUTHORCA025_GNM029/SOUTHORCA025_mesh_hgr.nc   mesh_hgr.nc
ln -s -v /store/njourd/SOUTHORCA025_GNM029/SOUTHORCA025_byte_mask.nc  mask.nc

####################################################################

for YEAR in $(seq $YEARi $YEARf)
do

# loop on files containing U,V
for file in ${INPUT_DIR}/${YEAR}/${CONFCASE}_y${YEAR}m??d??.5d_gridU.nc
do

  # input files for velocities :
  fileU=$file                                     # file containing uoce / vozocrtx
  fileV=`echo $file |sed -e "s/gridU/gridV/g"`    # file containing voce / vomecrty

  if [ $OPTION == '-ssh' ]; then
    fileT=`echo $file |sed -e "s/gridU/gridT/g"`  # file containing ssh / zos / sossheig
    OPTIONb="$OPTION $fileT"
  else
    OPTIONb="$OPTION"
  fi

  # output file containing the BSF :
  fileBSF=`echo $file |sed -e "s/gridU/psi/g"`

  ${CDFTOOL_DIR}/bin/cdfpsi -u $fileU -v $fileV $OPTIONb -o $fileBSF

  if [ -f $fileBSF ]; then
    echo "[oK]"
  else
    echo "~!@#%^&* ERROR: problem with creation of one of these files:"
    echo "                $fileBSF"
    echo "                >>>>>>>>>>>> stop !!"
    exit
  fi

done  # file

done  # YEAR
