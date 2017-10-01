#!/bin/bash

CONFIG="AMU12r"
BDY_DIR="$SHAREDELMER/input/nemo_${CONFIG}/BDY"
YEARi=1979
YEARf=1981

#----------------------------------------------------------------------------

for BDY in bdyT_tra bdyU_u3d bdyU_u2d bdyV_u3d bdyV_u2d bdyT_ice bdyT_ssh bdyU_ssh bdyV_ssh
do

for YEAR in $(seq $YEARi $YEARf)
do

  ncrcat ${BDY_DIR}/${BDY}_${YEAR}_*_${CONFIG}.nc ${BDY_DIR}/${BDY}_y${YEAR}_${CONFIG}.nc
  if [ -f ${BDY_DIR}/${BDY}_y${YEAR}_${CONFIG}.nc ]; then
    rm -f ${BDY_DIR}/${BDY}_${YEAR}_*_${CONFIG}.nc
    echo "${BDY_DIR}/${BDY}_y${YEAR}_${CONFIG}.nc  [oK]"
  else
    echo "~!@#%^&* ERROR: ${BDY_DIR}/${BDY}_y${YEAR}_${CONFIG}.nc HAS NOT BEEN CREATED   >>>>> STOP !!"
    exit
  fi

done

done
