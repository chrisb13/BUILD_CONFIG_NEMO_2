+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Nicolas Jourdain, IGE-CNRS, Feb. 2017
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Pre-processing tools to build a regional NEMO configuration laterally forced by a global ocean simulation.
Tested with NEMO-3.6 but should work with previous NEMO-3.x versions.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Before you start, make sure your system has the following:
- a fortran compiler
- netcdf and fortran-netcdf libraries.

To download NEMO, you need to register on http://www.nemo-ocean.eu/user/register

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Last updates:
- New version for general lon/lat grid (bi-linear interpolation)

Known caveats:
- Requires the same vertical discretization in the domain and for the datsaset used at BDYs and
  initial state (to be improved soon).

List of contributors: 
- Ute Hausmann

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#########################################################################################################
#########################################################################################################
## 0-- First of all, choose a configuration name (e.g. AMU12, WED12, PERIANT05) and clone this git repo :

	export CONFIG='AMU12.L75'  ## NB: to redo if you start a new session along the following steps !!

        cd $WORKDIR  ## can be your home, or any other directory...
        
        git clone https://github.com/nicojourdain/BUILD_CONFIG_NEMO_2.git BUILD_CONFIG_${CONFIG}
        # or : git clone git@github.com:nicojourdain/BUILD_CONFIG_NEMO_2.git BUILD_CONFIG_${CONFIG}
        #      (if you have a github account with your public key uploaded)


#########################################################################################################
#########################################################################################################
## 1-- Get and install NEMO and XIOS (required for the preprocessing to build the mesh_mask file)

        ## About NEMO/XIOS download and compilation: http://forge.ipsl.jussieu.fr/nemo/wiki/Users

        ## NB: on froggy (Ciment, Grenoble), you will need to run this command to enable svn functions:
        export ftp_proxy=http://www-cache.ujf-grenoble.fr:3128

        cd $WORKDIR  ## can be your home, or anywhere else
	mkdir models
	cd models

        ## Choose if you install XIOS-1 or XIOS-2 :
        svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0 XIOS  ## XIOS-1
        svn co -r 1011 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS     ## XIOS-2

        ## Install XIOS :
        cd XIOS
        chmod +x make_xios
        ./make_xios --avail      ## find your architecture if it exists (e.g. X64_ADA) 
                                 ## or create a new one in the arch directory.
        time ./make_xios --prod --full --arch X64_ADA --jobs 8 &> compil.log &
        tail -f compil.log       ## to follow the compilation
        ls bin/xios_server.exe   ## to check that it went ok
        cd ..

        ## Choose one of these NEMO releases (or any other branch and version) :
        svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk MY_NEMO         ## for the last trunk version
	svn co -r 6402 http://forge.ipsl.jussieu.fr/nemo/svn/trunk MY_NEMO ## for version 6402 of the trunk
        svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE MY_NEMO ## NEMO3.6 STABLE

        ## Compile for a new configuration, e.g. the AMU12 configuration (you can have several ones)
        cd MY_NEMO/NEMOGCM/CONFIG
        echo "$CONFIG OPA_SRC LIM_SRC_3" >> cfg.txt  ## here for OPA+LIM3 (see cfg.txt for other options)
        ./makenemo -h                ## find your architecture if it exists (e.g. X64_ADA)
                                     ## or create a new one in the MY_NEMO/NEMOGCM/ARCH directory.
        vi ../ARCH/arch-X64_ADA.fcm  ## adapt XIOS pathway accordingly with the previous installation.
        mkdir $CONFIG
        ## Choose fppkeys prior to compilation (see NEMO documentations or examples in existing configurations,
        ## e.g. MY_NEMO/NEMOGCM/CONFIG/AMM12/cpp_AMM12.fcm), e.g. : 
        echo " bld::tool::fppkeys key_bdy key_dynspg_ts key_lim3 key_ldfslp key_iomput key_mpp_mpi" > AMU12/cpp_${CONFIG}.fcm
        time ./makenemo -n $CONFIG -m X64_ADA -j 8 &> compile.log &
        tail -f compile.log            ## to follow the compilation
        ls -al ${CONFIG}/BLD/bin/nemo.exe  ## to check that everything went fine
        ## if you recompile, it is recommended to remove these directories: rm -rf ${CONFIG}/WORK ${CONFIG}/BLD
        cd ..

       ## Compile REBUILD_NEMO (used to recombine outputs at the end of the jobs):
       cd TOOLS
       ./maketools -h  ## use same architecture as for NEMO's compilation (e.g. X64_ADA)
       ./maketools -m X64_ADA -n REBUILD_NEMO
       ls -al REBUILD_NEMO/BLD/bin/  ## to check that everything went fine
       cd ..

#########################################################################################################
#########################################################################################################
## 2-- Prepare the pre-processing tools:

        ##### GSW TOOLBOX #####
        # First you may need the TEOS10 toolbox to convert from EOS80 to TEOS10.
        # To avoid issues with updates on the GSW-Fortran tools, I've cloned the 2016 GSW-Fortran
        # in this repository. In case you want to check for updates (not recommended due to important
        # updates), you can try with:    git clone https://github.com/TEOS-10/GSW-Fortran.git
        cd GSW-Fortran/test
        # edit makefile (in particular fill the FC and NETCDF_INCDIR variables)
        make
        ./gsw_check  ## to check (some have to work, maybe not all of them)
        cd ../..
        # The gsw_data_v3_0.nc file provided in this repository should be fine, but IF you start again 
        # from the original GSW git repository (i.e. from github.com/TEOS-10/GSW-Fortran.git), you will need to modify it
        # to avoid NaN in some places (otherwise, skip the next 4 lines):
        #   ./compile.sh remove_NaN_from_gsw_data_v3_0.f90 
        #   ./remove_NaN_from_gsw_data_v3_0
        #   ncks -x -v ocean_ref,ndepth_ref,deltaSA_ref,SA_ref,SAAR_ref GSW-Fortran/test/gsw_data_v3_0.nc gsw_data_v3_0.nc
        #   ncks -A gsw_data_v3_0_to_be_ncks-A.nc gsw_data_v3_0.nc 

        # Edit your own namelist, e.g. namelist_AMU12, and link namelist_pre to it:
        vi namelist_${CONFIG}  ## Adapt from provided namelist examples (check namelist comments).
        ln -s -v namelist_${CONFIG} namelist_pre  ## namelist_pre is the file read by all the scripts.
        # Edit compile_ALL.sh, adapt fortran compiler (and maybe netcdf path), then execute it:
        ./compile_ALL.sh
        # Adapt submit.sh to your platform (currently works for BATCH submission on occigen, ada, or 
        # for simple execution on Linux systems)

#########################################################################################################
#########################################################################################################
## 3-- Extract the coordinates files for the regional domain, build weights for interpolation, then
##     extract the bathymetry (and ice shelf draft if required):

        # fill the &general and &griddata section of the namelist

        ./submit.sh extract_coordinates 01  ## -> creates a "config_dir" directory (see namelist)
			                    ## -> creates the regional coordinates file 
                                            ##    coordinates_${CONFIG}.nc (stored in "config_dir")

        ./submit.sh build_linear_interp_weights 03  ## -> creates coeff_linear_${CONFIG}.nc (in "config_dir")

        # If you want to extract the bathymetry from a nemo simulation (e.g. eORCA12), while being 
        # consistent with the BDY bathymetry, do:
        ./submit.sh extract_bathy 01 30

        # Else, if you want to put RTOPO bathymetry within the domain while being consistent with the
        # BDY bathymetry, then fill the &bathy_special section of the namlelist and do:
        ./submit.sh extract_bathy_rtopo 04 30
        # Another option is to interpolate datasets from a stereographic grid (still with the 
        # &bathy_special section of the namlelist) :
        ./submit.sh extract_bathy_stereo 20 80

#########################################################################################################
#########################################################################################################
## 4- Create mesh_mask file.

        cd xxxxx ## TO BE IMPROVED
        ## NB: set jpni=jpnj=jpnij=0 in the &nammpp section of NEMO's namelist.
        ##     and use the same &namdom parameters as in the simulation used as BDYs.
        qsub run_nemo.sh    ## this will create mesh_mask_${CONFIG}.nc before crashing
        ./rebuild.sh mesh_mask
        mv mesh_mask.nc $WORKDIR/input/nemo_${CONFIG}/mesh_mask_${CONFIG}.nc ## directory defined as config_dir

#########################################################################################################
#########################################################################################################
## 5- Create the initial state (temperature and salinity)

        ./submit.sh extract_istate 01  ## -> should create the initial state for temperature and salinity 
                                       ##    e.g. dta_temp_AMU12.nc and dta_sal_AMU12.nc
                                       ##    (stored in directory defined as config_dir in the namelist)

#########################################################################################################
#########################################################################################################
## 5- Create the lateral boundary conditions (u, v, T, S, SSH, sea-ice thic, sea-ice frac)

	./submit.sh build_coordinates_bdy 01  ## -> creates the coordinate file for lateral boundaries
                                              ##    e.g. coordinates_bdy_AMU12.nc

        ./submit.sh extract_bdy_gridT 01 60   ## -> creates T,S bdy files and store them in the BDY folder
                                              ##    itself located in directory defined as config_dir

        ./submit.sh extract_bdy_icemod 01     ## -> creates ice bdy files and store them in the BDY folder

        ./submit.sh extract_bdy_ssh 01        ## -> creates SSH bdy files and store them in the BDY folder
                                              ##    NB: only bdyT_ssh are used by NEMO, but bdyU_ssh and 
                                              ##        bdyV_ssh are used by extract_bdy_psi.f90. 

        # Before extracting  baroclinic and barotropic velocities, you need to calculate the barotropic
        # stream function (BSF). To do so, install the CDFTOOLS: https://github.com/meom-group/CDFTOOLS
        # Then you can adapt compute_BSF.sh, and execute it:
        ./submit.sh 01 30 compute_BSF.sh      ## -> creates BSF files with same format as other inputs.

        # Then generate BDY files containing baroclinic velocities (simple interpolation):
        ./submit.sh extract_bdy_gridU 01 60   ## -> creates U3d bdy files and store them in the BDY folder
        ./submit.sh extract_bdy_gridV 01 60   ## -> creates V3d bdy files and store them in the BDY folder

        # Then generate BDY files containing barotropic velocities (from BSF to conserve transports):
        ./submit.sh extract_bdy_psi 01 20     ## -> creates U2d,V2d bdy files, store them in the BDY folder

        # To concatenate all the BDY files into yearly files, adapt concatenate_yearly_BDY.sh and execute: 
        ./concatenate_yearly_BDY.sh           ## -> concatenate to yearly files in the BDY folder
                                              ## NB: to do after extract_bdy_psi (ssh_U/V needed)

        # To prepare BDY files containing the tidal harmonics (only is you have tides at your BDYs) :
        ./submit.sh extract_bdy_tides 01 20   ## -> create bdytide files and store them in the BDY folder

#########################################################################################################
#########################################################################################################
## 6- Other files (SSS for restoring, runoff, chlorophyll)

	./submit.sh extract_SSS_restoring 01 30  ## -> creates SSS files and store them in a SSS folder
        ./concatenate_yearly_SSS.sh              ## Edit this file first.
                                                 ## -> concatenate the bdy files into yearly files

	./submit.sh extract_runoff_icebergs 01   ## -> creates iceberg runoff file
                                                 ##    (stored in config_dir defined in the namelist)

        ./submit.sh extract_chloro 01            ## -> creates chlorophyll file
                                                 ##    (stored in config_dir defined in the namelist)
    
#########################################################################################################
#########################################################################################################
## 7- Weights for the interpolation of atmospheric forcing

	# here it is assumed that you have NEMO downloaded and installed as explained in step #1

	cd $WORKDIR/models/MY_NEMO/NEMOGCM/TOOLS
	maketools -m X64_ADA -n WEIGHTS           ## Adapt for any machine different from ADA (try makenemo -h)
        cd WEIGHTS                                ## See nice README file in this directoy.
	cp -p nocsutil/namelist_example_bilin namelist_WEIGHTS_${CONFIG}_bilin ## Edit this namelist
	#NB: you may have to increase char_len in src/kinds_mod.f90 if you include long path for file names.

       ## Below is an example for namelist_WEIGHTS_${CONFIG}_bilin :
       ##      
       ##      &grid_inputs
       ##          input_file = 'drowned_precip_DFS5.2_y1993.nc'
       ##          nemo_file = 'coordinates_AMU12.nc'
       ##          datagrid_file = 'remap_data_grid.nc'
       ##          nemogrid_file = 'remap_nemo_grid.nc'
       ##          method = 'regular'
       ##          input_lon = 'lon'
       ##          input_lat = 'lat'
       ##          nemo_lon = 'glamt'
       ##          nemo_lat = 'gphit'
       ##          nemo_mask = 'none'
       ##          nemo_mask_value = 10
       ##          input_mask = 'none'
       ##          input_mask_value = 10
       ##      /
       ##      
       ##      &remap_inputs
       ##          num_maps = 1
       ##          grid1_file = 'remap_data_grid.nc'
       ##          grid2_file = 'remap_nemo_grid.nc'
       ##          interp_file1 = 'data_nemo_bilin.nc'
       ##          interp_file2 = 'nemo_data_bilin.nc'
       ##          map1_name = 'data to nemo bilin Mapping'
       ##          map2_name = 'nemo to data bilin Mapping'
       ##          map_method = 'bilinear'
       ##          normalize_opt = 'frac'
       ##          output_opt = 'scrip'
       ##          restrict_type = 'latitude'
       ##          num_srch_bins = 90
       ##          luse_grid1_area = .false.
       ##          luse_grid2_area = .false.
       ##      /
       ##      
       ##      &interp_inputs
       ##          input_file = "drowned_precip_DFS5.2_y1993.nc"
       ##          interp_file = "data_nemo_bilin.nc"
       ##          input_name = "snow"
       ##          input_start = 1,1,1,1
       ##          input_stride = 1,1,1,1
       ##          input_stop = 0,0,0,1
       ##          input_vars = 'initial_time0_hours'
       ##      /
       ##      
       ##      &interp_outputs
       ##          output_file = "snow_orca.nc"
       ##          output_mode = "create"
       ##          output_dims = 'x', 'y', 'time_counter'
       ##          output_scaling = "snow|1.0", "time_counter|86400.0"
       ##          output_name = 'snow'
       ##          output_lon = 'x'
       ##          output_lat = 'y'
       ##          output_vars = 'time_counter'
       ##          output_attributes = 'time_counter|units|seconds since 1995-00-00 00:00:00',
       ##                              'time_counter|calendar|noleap',
       ##                              'snow|units|mm/s'
       ##      /
       ##      

       ## Then, execute these things:
	./scripgrid.exe   ## namelist_${CONFIG}_bilin
	./scrip.exe       ## namelist_${CONFIG}_bilin
	./scripshape.exe  ## namelist_${CONFIG}_bilin

       ## Then create another namelist:  namelist_${CONFIG}_bicub and just change these parameters:
       ##      interp_file1 = 'data_nemo_bicubic.nc'
       ##      interp_file2 = 'nemo_data_bicubic.nc'
       ##      map1_name = 'data to nemo bicubic Mapping'
       ##      map2_name = 'nemo to data bicubic Mapping'
       ##      map_method = 'bicubic'
       ##      interp_file = "data_nemo_bicubic.nc"
       ##      interp_file = 'data_nemo_bicubic.nc'
       ##      output_file = 'weights_bicubic.nc'
       ## Then execute :
	./scrip.exe       ## namelist_${CONFIG}_bicub
	./scripshape.exe  ## namelist_${CONFIG}_bicub
