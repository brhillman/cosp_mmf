! ------------------------------------------------------------------------------
!
! Driver to run COSP on MMF output
!
! Authors: Benjamin R. Hillman, Roger Marchand
!  (BRH adapted from original driver written by Roger Marchand)
!
! History:
!  Apr 2013: Roger Marchand modified driver and COSP to work with new MMF V5, 
!            which can use an adaptive vertical grid
!  Aug 2015: Benjamin Hillman simplified a bit; removed options to modify CRM
!            inputs. These are now handled externally to improve transparency
!            in workflow, making this driver a more generic piece of code that
!            is only designed to read in MMF output and run COSP.
! 
! ------------------------------------------------------------------------------
program cosp_mmf
   use cosp_mmf_utils
   use mod_cosp_constants
   use mod_cosp_io
   use mod_cosp_types
   use mod_cosp

   implicit none

   ! MMF input variables
   type(mmf_input_type) :: mmf_input

   ! COSP input namelist variables
   character(len = 64) :: cosp_input_nl = 'cosp_input_nl.txt'
   character(len = 64) :: cosp_output_nl = 'cosp_output_nl.txt'
   character(len = 64) :: cmor_nl
   integer :: overlap
   integer :: isccp_topheight, isccp_topheight_direction
   integer :: npoints, npoints_it, ncolumns, nlevels, nlr
   logical :: use_vgrid, csat_vgrid
   real :: radar_freq, k2
   integer :: surface_radar, use_mie_tables, use_gas_abs
   integer :: do_ray, melt_lay
   integer :: nprmts_max_hydro, naero, nprmts_max_aero
   integer :: lidar_ice_type
   logical :: use_precipitation_fluxes, use_reff
   integer :: platform, satellite, instrument, nchannels
   integer, dimension(RTTOV_MAX_CHANNELS) :: channels
   real, dimension(RTTOV_MAX_CHANNELS) :: surfem
   real :: zenang
   real :: co2, ch4, n2o, co

   ! COSP gridbox variables
   real :: emsfc_lw = 0.99

   ! COSP point outputs
   type(cosp_config) :: cfg ! options
   type(cosp_gridbox) :: gbx ! gridbox inputs
   type(cosp_subgrid) :: sgx ! subgrid outputs
   type(cosp_sgradar) :: sgradar ! output from radar simulator
   type(cosp_sglidar) :: sglidar ! output from lidar simulator
   type(cosp_isccp) :: isccp ! output from ISCCP simulator
   type(cosp_modis) :: modis ! output from MODIS simulator
   type(cosp_misr) :: misr ! output from MISR simulator
   type(cosp_vgrid) :: vgrid ! information on vertical grid of stats
   type(cosp_radarstats) :: stradar ! summary statistics from radar simulator
   type(cosp_lidarstats) :: stlidar ! summary statistics from lidar simulator

   ! COSP global outputs
   type(cosp_gridbox) :: gbx_globe ! gridbox inputs
   type(cosp_subgrid) :: sgx_globe ! subgrid outputs
   type(cosp_sgradar) :: sgradar_globe ! output from radar simulator
   type(cosp_sglidar) :: sglidar_globe ! output from lidar simulator
   type(cosp_isccp) :: isccp_globe ! output from ISCCP simulator
   type(cosp_modis) :: modis_globe ! output from MODIS simulator
   type(cosp_misr) :: misr_globe ! output from MISR simulator
   type(cosp_vgrid) :: vgrid_globe ! information on vertical grid of stats
   type(cosp_radarstats) :: stradar_globe ! summary statistics from radar simulator
   type(cosp_lidarstats) :: stlidar_globe ! summary statistics from lidar simulator

   ! added to get cmor output from development version of COSP (Sept 2013 R82)
   type(cosp_rttov) :: rttov, rttov_globe ! output from RTTOV
   type(var1d) :: v1d(N1D) ! structures needed by output routines for 1D variables
   type(var2d) :: v2d(N2D) ! structures needed by output routines for 2D variables
   type(var3d) :: v3d(N3D) ! structures needed by output routines for 3D variables
   integer :: grid_id, latvar_id, lonvar_id, lon_axid, lat_axid, time_axid, height_axid, height_mlev_axid, column_axid, sza_axid, &
               temp_axid, channel_axid, dbze_axid, sratio_axid, MISR_CTH_axid, tau_axid, pressure2_axid

   ! local variables
   character(len = 512) :: inputfile, arg
   integer :: idx
   integer :: nlev, ntime
   integer :: ilon, ilat, itime, ilev, ihydro, ipoint, icol
   integer :: npoints_globe, ncolumns_globe
   logical :: verbose = .false.

   ! COSP input namelist
   namelist/cosp_input/cmor_nl, overlap, isccp_topheight, &
         isccp_topheight_direction, &
         npoints, npoints_it, ncolumns, nlevels, use_vgrid, nlr, &
         csat_vgrid, &
         radar_freq, surface_radar, use_mie_tables, &
         use_gas_abs, do_ray, melt_lay, k2, &
         nprmts_max_hydro, naero, nprmts_max_aero, &
         lidar_ice_type, use_precipitation_fluxes, use_reff, &
         platform, satellite, instrument, nchannels, &
         channels, surfem, zenang, co2, ch4, n2o, co

   ! get filename from command line arguments
   if (iargc() == 0) then
      print *, 'usage: cosp_driver <input file>'
      stop
   else
      call getarg(1, inputfile)
   end if

   ! read data from input file
   if (verbose) print *, 'Read MMF input file...'
   call read_mmf_input(inputfile, mmf_input)      

   ! read COSP input namelist
   open(10, file = cosp_input_nl, status = 'old')
   read(10, nml = cosp_input)
   close(10)

   ! read COSP output namelist
   call read_cosp_output_nl(cosp_output_nl, cfg)

   ! set ncolumns and npoints (override cosp_input_nl)
   ncolumns = 1
   npoints = mmf_input%nx * mmf_input%ny

   ! set nlevels as the maximum number of levels that will used
   ! (overrides cosp_input_nl.txt)
   if(mmf_input%azg_input) then
      nlevels = mmf_input%nz
   else
      nlevels = mmf_input%nlev
   endif

   ! loop over time samples
   ! NOTE: this has not been tested with more than one time sample!!!
   do itime = 1, mmf_input%ntime
      ! initialize output structures for global output
      if (verbose) print *, 'Initialize COSP global types...'
      npoints_globe = mmf_input%nlon * mmf_input%nlat
      ncolumns_globe = npoints
      call construct_cosp_gridbox( &
         mmf_input%time(itime), mmf_input%time_bnds(itime, :), &
         radar_freq, surface_radar, &
         use_mie_tables, use_gas_abs, &
         do_ray, melt_lay, k2, &
         npoints_globe, nlevels, ncolumns_globe, N_HYDRO, &
         nprmts_max_hydro, Naero, &
         nprmts_max_aero, npoints_it, &
         lidar_ice_type, isccp_topheight, &
         isccp_topheight_direction, &
         overlap, emsfc_lw, &
         use_precipitation_fluxes, use_reff, &
         platform, satellite, instrument, nchannels, zenang, &
         channels(1:nchannels), surfem(1:nchannels), &
         co2, ch4, n2o, co, &
         gbx_globe, load_LUT = .true. &
      )
      call construct_cosp_vgrid( &
         gbx_globe, Nlr, use_vgrid, csat_vgrid, &
         vgrid_globe &
      )
      call construct_cosp_subgrid( &
         npoints_globe, ncolumns_globe, nlevels, &
         sgx_globe &
      )
      call construct_cosp_sgradar( &
         cfg, npoints_globe, ncolumns_globe, nlevels, N_HYDRO, &
         sgradar_globe &
      )
      call construct_cosp_radarstats( &         
         cfg, npoints_globe, ncolumns_globe, vgrid_globe%nlvgrid, &
         N_HYDRO, &
         stradar_globe &
      )
      call construct_cosp_sglidar( &         
         cfg, npoints_globe, ncolumns_globe, nlevels, &
         N_HYDRO, PARASOL_NREFL, &
         sglidar_globe &
      )
      call construct_cosp_lidarstats( &         
         cfg, npoints_globe, ncolumns_globe, vgrid_globe%nlvgrid, &
         N_HYDRO, PARASOL_NREFL, &
         stlidar_globe &
      )
      call construct_cosp_isccp( &
         cfg, npoints_globe, ncolumns_globe, nlevels, &
         isccp_globe &
      )
      call construct_cosp_modis(cfg, npoints_globe, modis_globe)
      call construct_cosp_misr(cfg, npoints_globe, misr_globe)
      call construct_cosp_rttov(cfg, npoints, nchannels, rttov_globe)

      ! loop over lat/lon and run COSP separately for each gridbox
      if (verbose) print *, 'Loop over lon/lat; run cosp...'
      do ilon = 1, mmf_input%nlon
      do ilat = 1, mmf_input%nlat

         ! map 2D lat/lon index to 1D gridbox index 
         ipoint = (ilon + (ilat - 1) * mmf_input%nlon)

         ! if using adaptive grid we reset nlevels used by COSP equal to 
         ! number of layers used in this grid-cell
         if(mmf_input%azg_input) then
            nlevels = mmf_input%nzm_used(ilon, ilat, itime)
         endif

         ! initialize local gridbox
         if (ipoint == 1 .or. mmf_input%azg_input) then
            if(ipoint > 1) then  
               ! this is not first pass of loop ...
               ! need to clear previous local gridbox
               call free_cosp_gridbox(gbx, save_LUT = .false.)
               call free_cosp_subgrid(sgx)
               call free_cosp_sgradar(sgradar)
               call free_cosp_radarstats(stradar)
               call free_cosp_sglidar(sglidar)
               call free_cosp_lidarstats(stlidar)
               call free_cosp_isccp(isccp)
               call free_cosp_misr(misr)
               call free_cosp_modis(modis)
               call free_cosp_rttov(rttov)
               call free_cosp_vgrid(vgrid)
            endif

            call construct_cosp_gridbox( &
               mmf_input%time(itime), mmf_input%time_bnds(itime, :), &
               radar_freq, surface_radar, &
               use_mie_tables, use_gas_abs, &
               do_ray, melt_lay, k2, &
               npoints, nlevels, ncolumns, N_HYDRO, &
               nprmts_max_hydro, Naero, &
               nprmts_max_aero, npoints_it, &
               lidar_ice_type, isccp_topheight, &
               isccp_topheight_direction, &
               overlap, emsfc_lw, &
               use_precipitation_fluxes, use_reff, &
               platform, satellite, instrument, nchannels, zenang, &
               channels(1:nchannels), surfem(1:nchannels), &
               co2, ch4, n2o, co, &
               gbx, load_LUT = .false. &
            )

            ! NOTE: we do not load LUT when creating the local grid box ...
            !       simply copy from the global structure
            gbx%hp%Z_scale_flag = gbx_globe%hp%Z_scale_flag
            gbx%hp%Ze_scaled = gbx_globe%hp%Ze_scaled
            gbx%hp%Zr_scaled = gbx_globe%hp%Zr_scaled
            gbx%hp%kr_scaled = gbx_globe%hp%kr_scaled  

            ! initialize outputs for local gridbox
            call construct_cosp_vgrid(gbx, Nlr, use_vgrid, csat_vgrid, vgrid)
            call construct_cosp_subgrid(npoints, ncolumns, nlevels, sgx)
            call construct_cosp_sgradar(cfg, npoints, ncolumns, nlevels, &
                                        N_HYDRO, sgradar)
            call construct_cosp_radarstats(cfg, npoints, ncolumns, &
                                           vgrid%nlvgrid, N_HYDRO, &   
                                           stradar)
            call construct_cosp_sglidar(cfg, npoints, ncolumns, nlevels, &   
                                        N_HYDRO, PARASOL_NREFL, sglidar)
            call construct_cosp_lidarstats(cfg, npoints, ncolumns, &
                                           vgrid%nlvgrid, &   
                                           N_HYDRO, PARASOL_NREFL, stlidar) 
            call construct_cosp_isccp(cfg, npoints, ncolumns, nlevels, isccp)
            call construct_cosp_modis(cfg, npoints, modis)
            call construct_cosp_misr(cfg, npoints, misr)
            call construct_cosp_rttov(cfg, npoints, nchannels, rttov)
         end if ! construct local COSP gbx, if needed

         ! populate COSP input gridbox with MMF outputs
         call populate_gridbox(mmf_input, gbx, ilon, ilat, itime)

         ! run COSP
         call cosp(overlap, ncolumns, cfg, vgrid, gbx, sgx, &
                   sgradar, sglidar, isccp, misr, modis, stradar, stlidar)

         ! save COSP point outputs in global structure
         call save_output( &
             cfg, &
             gbx, sgx, sgradar, sglidar, &
             isccp, misr, modis, stradar, stlidar, &
             gbx_globe, sgx_globe, sgradar_globe, sglidar_globe, &
             isccp_globe, misr_globe, modis_globe, &
             stradar_globe, stlidar_globe, &
             loc = ipoint &
         )

         ! save current LUT
         gbx_globe%hp%Z_scale_flag = gbx%hp%Z_scale_flag
         gbx_globe%hp%Ze_scaled = gbx%hp%Ze_scaled
         gbx_globe%hp%Zr_scaled = gbx%hp%Zr_scaled
         gbx_globe%hp%kr_scaled = gbx%hp%kr_scaled  
      end do ! loop over lat
      end do ! loop over lon

      ! write COSP outputs
      if (cfg%Lwrite_output) then
         if (verbose) print *, 'Write COSP global outputs...'
         call nc_cmor_init( &
            cmor_nl, 'replace', cfg, &
            vgrid_globe, gbx_globe, sgx_globe, sgradar_globe, sglidar_globe, &
            isccp_globe, misr_globe, modis_globe, &
            rttov_globe, stradar_globe, stlidar_globe, &
            2, mmf_input%nlon, mmf_input%nlat, N1D, N2D, N3D, &
            lon_axid, lat_axid, time_axid, height_axid, &
            height_mlev_axid, grid_id, lonvar_id, latvar_id, &
            column_axid, sza_axid, temp_axid, channel_axid, dbze_axid, &
            sratio_axid, MISR_CTH_axid, tau_axid, pressure2_axid, &
            v1d, v2d, v3d &
         )
         call nc_cmor_write_2d( &
                  mmf_input%time(itime), mmf_input%time_bnds(itime, :), &
                  2, mmf_input%nlon, mmf_input%nlat, &
                  N1D, N2D, N3D, v1d, v2d, v3d &
         )
         call nc_cmor_close()
      endif

      ! deallocate COSP point outputs
      if (verbose) print *, 'Deallocate COSP gridbox types...'
      call free_cosp_gridbox(gbx, save_LUT = .false.)
      call free_cosp_subgrid(sgx)
      call free_cosp_sgradar(sgradar)
      call free_cosp_radarstats(stradar)
      call free_cosp_sglidar(sglidar)
      call free_cosp_lidarstats(stlidar)
      call free_cosp_isccp(isccp)
      call free_cosp_misr(misr)
      call free_cosp_modis(modis)
      call free_cosp_rttov(rttov)
      call free_cosp_vgrid(vgrid)

      ! deallocate COSP global outputs
      if (verbose) print *, 'Deallocate COSP global types...'
      call free_cosp_gridbox(gbx_globe, save_LUT = .true.)
      call free_cosp_subgrid(sgx_globe)
      call free_cosp_sgradar(sgradar_globe)
      call free_cosp_radarstats(stradar_globe)
      call free_cosp_sglidar(sglidar_globe)
      call free_cosp_lidarstats(stlidar_globe)
      call free_cosp_isccp(isccp_globe)
      call free_cosp_misr(misr_globe)
      call free_cosp_modis(modis_globe)
      call free_cosp_rttov(rttov_globe) ! R82
      call free_cosp_vgrid(vgrid_globe)
   end do ! loop over time indices

   ! deallocate MMF fields
   if (verbose) print *, 'Deallocate MMF fields...'
   call deallocate_mmf_input(mmf_input)

   if(verbose) print *, 'Done.'
end program cosp_mmf
