! ------------------------------------------------------------------------------
!
! Module containing derived types and subroutines to work with MMF output
!
! Authors: Benjamin R. Hillman, Roger Marchand
!
! History:
!  Original version: BRH adapted from original driver written by Roger Marchand
!  April 2013: Roger Marchand modified to work with output from MMF V5 with
!              adaptive vertical grid
!  May 2015: BRH modified significantly
!            New routines to mimic GCM assumptions using MMF CRM fields
!  August 2015: BRH simplified
!            Removed routines to mimic GCM assumptions; these now handled
!            externally for greater transparency in workflow
!
! ------------------------------------------------------------------------------
module cosp_mmf_utils
   use netcdf
   use netcdf_utils
   use mod_cosp_constants
   use mod_cosp_types
   use mod_cosp_modis_simulator
   use mod_cosp_io
   use shr_kind_mod, only: r8=>shr_kind_r8
   use pkg_cldoptics, only: cldefr, cldems

   implicit none
   private
   public mmf_input_type, allocate_mmf_input, deallocate_mmf_input, &
          read_mmf_dimensions, read_mmf_input, populate_gridbox, save_output

   ! ---------------------------------------------------------------------------
   ! derived type to hold fields from MMF output files
   ! ---------------------------------------------------------------------------
   type mmf_input_type

      ! dimension sizes (BRH added 2013/04/30)
      integer :: ntime
      integer :: nlat
      integer :: nlon
      integer :: nlev
      integer :: nx
      integer :: ny
      integer :: nz

      ! flag to hold initialized state of structure
      logical :: initialized = .false.

      ! flag for adaptive grid version (BRH added 2015/05/19)
      logical :: azg_input = .false.

      ! GCM inputs in MMF3/5, but CRM inputs in MMF5_azg (adaptive grid version)
      real,pointer :: pmid(:,:,:,:)  ! pressure at midpoints
      real,pointer :: pint(:,:,:,:)  ! pressure at midpoints
      real,pointer :: zmid(:,:,:,:)    ! height in meters
      real,pointer :: zint(:,:,:,:)    ! height in meters
      real,pointer :: t(:,:,:,:)    ! temperature in kelvin
      real,pointer :: sh(:,:,:,:) ! specific humidity in kg/kg
      real,pointer :: rh(:,:,:,:) ! relative humidity

      ! GCM inputs
      real,pointer :: psfc(:,:,:) ! surface pressure
      real,pointer :: tsfc(:,:,:) ! surface temperature
      real,pointer :: sunlit(:,:,:)    ! solar insolation
      double precision,pointer :: lev(:) ! nominal pressure levels (hybrid?)
      double precision,pointer :: lat(:) ! latitude
      double precision,pointer :: lon(:) ! longitude
      double precision,pointer :: time(:) ! time
      double precision,pointer :: time_bnds(:,:) ! time bounds

      ! CRM inputs
      double precision,pointer :: crm_z(:) ! nominal CRM pressure levels
      real,pointer :: nzm_used(:,:,:)   ! Number of layer (mid-points) used of nz = nz_max possible ... added for MMF5_azg
      real,pointer :: qc(:,:,:,:,:)   ! cloud non-precipitating water
      real,pointer :: qi(:,:,:,:,:)   ! cloud non-precipitating ice
      real,pointer :: qpc(:,:,:,:,:)  ! cloud precipitating water
      real,pointer :: qpi(:,:,:,:,:)  ! cloud precipitating ice
      real,pointer :: qps(:,:,:,:,:)  ! cloud precipitating ice
      real,pointer :: qpg(:,:,:,:,:)   ! graupel
      real,pointer :: nc(:,:,:,:,:)   ! cloud water number concentration
      real,pointer :: ni(:,:,:,:,:)   ! cloud ice number concentration
      real,pointer :: nr(:,:,:,:,:)   ! cloud precipitating water (rain) number concentration
      real,pointer :: ns(:,:,:,:,:)   ! cloud precipitating ice (snow) number concentration
      real,pointer :: ng(:,:,:,:,:)   ! cloud precipitating ice (graupel) number concentration
   end type mmf_input_type

   logical :: verbose = .false. ! BRH added global verbose flag 2013/04/30

contains

   ! ---------------------------------------------------------------------------
   ! subroutine to allocate memory for MMF fields
   ! ---------------------------------------------------------------------------
   subroutine allocate_mmf_input(mmf_input)
      type(mmf_input_type),intent(inout) :: mmf_input

      ! local variables
      integer :: ntime,nlat,nlon,nlev,nx,ny,nz

      ! dimension sizes (BRH 2013/04/30)
      ntime = mmf_input%ntime
      nlat = mmf_input%nlat
      nlon = mmf_input%nlon
      nlev = mmf_input%nlev
      nx = mmf_input%nx
      ny = mmf_input%ny
      nz = mmf_input%nz

      ! check dimensions
      if (ntime == -1 .or. nlat == -1 .or. nlon == -1 .or. nlev == -1 .or. nx == -1 .or. ny == -1 .or. nz == -1) then 
         print *, 'One or more dimensions not found; please call read_mmf_dimensions first.'
         stop
      end if

      ! vertical pressure/height grid may be on GCM or CRM grid
      if (mmf_input%azg_input) then
         allocate(mmf_input%lev(nz))
         allocate(mmf_input%pmid(nlon,nlat,nz,ntime))
         allocate(mmf_input%pint(nlon,nlat,nz,ntime))
         allocate(mmf_input%zmid(nlon,nlat,nz,ntime))
         allocate(mmf_input%zint(nlon,nlat,nz,ntime))
         allocate(mmf_input%t(nlon,nlat,nz,ntime))
         allocate(mmf_input%sh(nlon,nlat,nz,ntime))
         allocate(mmf_input%rh(nlon,nlat,nz,ntime))
      else
         allocate(mmf_input%lev(nlev))
         allocate(mmf_input%pmid(nlon,nlat,nlev,ntime))
         allocate(mmf_input%pint(nlon,nlat,nlev,ntime))
         allocate(mmf_input%zmid(nlon,nlat,nlev,ntime))
         allocate(mmf_input%zint(nlon,nlat,nlev,ntime))
         allocate(mmf_input%t(nlon,nlat,nlev,ntime))
         allocate(mmf_input%sh(nlon,nlat,nlev,ntime))
         allocate(mmf_input%rh(nlon,nlat,nlev,ntime))
      end if

      ! allocate memory for GCM fields
      allocate(mmf_input%sunlit(nlon,nlat,ntime))
      allocate(mmf_input%psfc(nlon,nlat,ntime))
      allocate(mmf_input%tsfc(nlon,nlat,ntime))
      allocate(mmf_input%lat(nlat))
      allocate(mmf_input%lon(nlon))
      allocate(mmf_input%time(ntime))
      allocate(mmf_input%time_bnds(ntime,2))

      ! initialize GCM fields
      mmf_input%sunlit = 0
      mmf_input%psfc = 0
      mmf_input%tsfc = 0
      mmf_input%pmid = 0
      mmf_input%pint = 0
      mmf_input%zmid = 0
      mmf_input%zint = 0
      mmf_input%t = 0
      mmf_input%sh = 0
      mmf_input%rh = 0
      mmf_input%lev = 0
      mmf_input%lat = 0
      mmf_input%lon = 0
      mmf_input%time = 0

      ! allocate memory for CRM fields
      allocate(mmf_input%crm_z(nlev))
      allocate(mmf_input%nzm_used(nlon,nlat,ntime)) ! Roj added for MMF5_azg
      allocate(mmf_input%qc(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%qi(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%qpc(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%qpi(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%qps(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%qpg(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%nc(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%ni(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%nr(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%ns(nlon,nlat,nx*ny,nlev,ntime))
      allocate(mmf_input%ng(nlon,nlat,nx*ny,nlev,ntime))

      ! initialize CRM fields
      mmf_input%crm_z = 0
      mmf_input%nzm_used = 0
      mmf_input%qc = 0
      mmf_input%qi = 0
      mmf_input%qpc = 0
      mmf_input%qpg = 0
      mmf_input%qpi = 0
      mmf_input%qps = 0
      mmf_input%nc = 0
      mmf_input%ni = 0
      mmf_input%nr = 0
      mmf_input%ns = 0
      mmf_input%ng = 0

      mmf_input%initialized = .true.
   end subroutine allocate_mmf_input


   ! ---------------------------------------------------------------------------
   ! subroutine to deallocate memory for MMF fields
   ! ---------------------------------------------------------------------------
   subroutine deallocate_mmf_input(mmf_input)
      type(mmf_input_type),intent(inout) :: mmf_input

      ! deallocate GCM fields
      deallocate(mmf_input%sunlit)
      deallocate(mmf_input%psfc)
      deallocate(mmf_input%tsfc)
      deallocate(mmf_input%pmid)
      deallocate(mmf_input%pint)
      deallocate(mmf_input%zmid)
      deallocate(mmf_input%zint)
      deallocate(mmf_input%t)
      deallocate(mmf_input%sh)
      deallocate(mmf_input%rh)
      deallocate(mmf_input%lev)
      deallocate(mmf_input%lat)
      deallocate(mmf_input%lon)
      deallocate(mmf_input%time)
      deallocate(mmf_input%time_bnds)

      ! deallocate CRM fields
      deallocate(mmf_input%crm_z)    ! added by Ben
      deallocate(mmf_input%nzm_used) ! added by Roj
      deallocate(mmf_input%qc)
      deallocate(mmf_input%qi)
      deallocate(mmf_input%qpc)
      deallocate(mmf_input%qpi)
      deallocate(mmf_input%qps)
      deallocate(mmf_input%qpg)
      deallocate(mmf_input%nc)
      deallocate(mmf_input%ni)
      deallocate(mmf_input%nr)
      deallocate(mmf_input%ns)
      deallocate(mmf_input%ng)

      mmf_input%initialized = .false.
   end subroutine deallocate_mmf_input


   ! ---------------------------------------------------------------------------
   ! subroutine to read MMF dimensions from file
   ! ---------------------------------------------------------------------------
   subroutine read_mmf_dimensions(inputfile,mmf_input)
      character(len=*),intent(in) :: inputfile
      type(mmf_input_type),intent(inout) :: mmf_input
      integer :: ncid

      ! open file
      call nc_check(nf90_open(inputfile,nf90_nowrite,ncid))

      ! get large-scale grid dimensions
      call nc_get_dim(ncid,'lon',mmf_input%nlon)
      call nc_get_dim(ncid,'lat',mmf_input%nlat)
      call nc_get_dim(ncid,'lev',mmf_input%nlev)
      call nc_get_dim(ncid,'time',mmf_input%ntime)

      ! find CRM grid dimensions
      mmf_input%azg_input = .false.  ! Initially assume this is not output
                                     ! from the adaptive grid version of MMF5 
      call nc_get_dim(ncid,'crm_x',mmf_input%nx)   
                                     ! In original MMF output files the crm x 
                                     ! variable dimension is "crm_x", 
                                     ! but changed to "crm_nx" in MMF5
      if(mmf_input%nx == -1) then    ! did NOT find crm_x
         if (verbose) print *,'   crm_x not found ... will look for crm_nx'
         call nc_get_dim(ncid,'crm_nx',mmf_input%nx)
         call nc_get_dim(ncid,'crm_ny',mmf_input%ny)     

         ! Is this an adaptive grid version of MMF5 ? 
         ! ... look for "crm_nz_max", which indicates adaptive grid version of MMF
         call nc_get_dim(ncid,'crm_nzm_max',mmf_input%nz)
         if(mmf_input%nz == -1) then ! did NOT find crm_nzm_max
            if (verbose) print *, '   crm_nzm_max not found ... this is output from MMF5 without an adaptive grid' 
            ! No. then this is an MMF5 version w/out adaptive grid so "crm_nz" should be defined
            call nc_get_dim(ncid,'crm_nz',mmf_input%nz)
         else
            if (verbose) print *, '   crm_nzm_max found ... this is output from MMF5 with adaptive grid'
            mmf_input%azg_input = .true.
         endif
      else
         ! if we did find "crm_x" we expect "crm_y" and "crm_z" ...
         call nc_get_dim(ncid,'crm_y',mmf_input%ny)
         call nc_get_dim(ncid,'crm_z',mmf_input%nz)
      endif

      ! make sure we found all dimensions
      if (mmf_input%nz == -1 .or. mmf_input%nx == -1 .or. mmf_input%ny == -1 &
            .or. mmf_input%nlon == -1 .or. mmf_input%nlat == -1 &
            .or. mmf_input%nlev == -1 .or. mmf_input%ntime == -1 &
      ) then 
         print *,'   Unable to find all necessary dimensions lengths.'
         stop
      endif

      ! close file
      call nc_check(nf90_close(ncid))
   end subroutine read_mmf_dimensions


   ! ---------------------------------------------------------------------------
   ! subroutine to read MMF fields from file
   ! ---------------------------------------------------------------------------
   subroutine read_mmf_input(inputfile,mmf_input)
      character(len = *), intent(in) :: inputfile
      type(mmf_input_type), intent(out) :: mmf_input

      ! local variables
      integer :: nlon, nlat, nlev, nx, ny, nz, ntime
                        ! Note: these are stored as part of mmf_input structure 
                        ! and don't need to be passed as subroutine parameters
                        ! anymore (BRH 2015/05/19)
      integer :: i, ilon, ilat, itime
      integer :: ix, iy, iz
      real,allocatable :: phis(:,:,:), pint(:,:,:,:), dp(:,:,:,:)
      real,allocatable :: wg_ratio(:,:,:,:,:)
      real,allocatable :: v6d(:,:,:,:,:,:)
      character*64 :: attval
      integer :: ncid, varid
      logical :: gcm_flipped = .false.
      logical :: crm_flipped = .false.

      ! make sure structures are allocated
      if (.not. mmf_input%initialized) then
         call read_mmf_dimensions(inputfile, mmf_input)
         call allocate_mmf_input(mmf_input)
      end if

      ! open file
      call nc_check(nf90_open(inputfile, nf90_nowrite, ncid))

      ! copy dimension sizes from mmf_input structure for convenience
      ntime = mmf_input%ntime
      nlon = mmf_input%nlon
      nlat = mmf_input%nlat
      nlev = mmf_input%nlev
      nx = mmf_input%nx
      ny = mmf_input%ny
      nz = mmf_input%nz

      ! coordinate variables
      call nc_get_data(ncid, 'crm_z', d1d=mmf_input%crm_z)
      call nc_get_data(ncid, 'lev', d1d = mmf_input%lev)
      call nc_get_data(ncid, 'lat', d1d = mmf_input%lat)
      call nc_get_data(ncid, 'lon', d1d = mmf_input%lon)
      call nc_get_data(ncid, 'time', d1d = mmf_input%time)

      ! We have to add to the time coordinate because the MMF
      ! "days since" time does not match up with the calendar
      ! date for some reason...why is this?
      ! BRH moved this to subroutine and fixed 2015-05-22...
      ! previous versions failed to fix the time coordinate,
      ! so dates in COSP output files did not match the dates
      ! of the MMF input files. They should be consistent now though.
      mmf_input%time = mmf_input%time + 1

      ! PNNL MMF files have weird time bounds for some reason,
      ! and CMOR changes the time coordinate to be the midpoint
      ! between the two bounds, so I define my own bounds here
      ! to straddle the time coordinate.
      mmf_input%time_bnds(:,1) = mmf_input%time(:) - 0.125 * 0.5
      mmf_input%time_bnds(:,2) = mmf_input%time(:) + 0.125 * 0.5

      ! surface fields
      call nc_get_data(ncid, 'PS', v3d = mmf_input%psfc)
      call nc_get_data(ncid, 'TS', v3d = mmf_input%tsfc)
      call nc_get_data(ncid, 'SOLIN', v3d = mmf_input%sunlit)
      where (mmf_input%sunlit > 0)
         mmf_input%sunlit = 1
      elsewhere
         mmf_input%sunlit = 0
      endwhere

      ! level fields
      if (mmf_input%azg_input) then  ! is this adaptive grid output?
                                     ! if so should have variable CRM_NZM
         if (verbose) print *, 'MMF5 AZG input file; use CRM p, dp, z' 

         call nc_get_data(ncid, 'CRM_NZM', v3d = mmf_input%nzm_used)
         call nc_get_data(ncid, 'CRM_PM',  v4d = mmf_input%pmid)
         call nc_get_data(ncid, 'CRM_PI',  v4d = mmf_input%pint(:,:,:nz,:))
         call nc_get_data(ncid, 'CRM_ZM',  v4d = mmf_input%zmid)

         ! add surface height to height above surface
         allocate(phis(nlon,nlat,ntime))
         call nc_get_data(ncid, 'PHIS', v3d = phis)
         do iz = 1, nz
            mmf_input%zmid(:,:,iz,:) = mmf_input%zmid(:,:,iz,:) &
                + phis(:,:,:)/9.8
         end do
         deallocate(phis)

      else ! use pressure/height grid from GCM 
         if (verbose) print *, 'MMF3 input file; use GCM p, dp, z' 

         ! pressure
         allocate(dp(nlon,nlat,nlev,ntime))
         call nc_get_data(ncid,'PRES',v4d = mmf_input%pmid)
         call nc_get_data(ncid,'DPRES',v4d = dp)
         mmf_input%pint = mmf_input%pmid + 0.5 * dp
         deallocate(dp)

         ! height
         call nc_get_data(ncid, 'Z3', v4d = mmf_input%zmid)
         mmf_input%zint(:,:,1,:) = mmf_input%zmid(:,:,1,:) - 0.5 * abs( &
            mmf_input%zmid(:,:,2,:) - mmf_input%zmid(:,:,1,:) &
         )
         mmf_input%zint(:,:,2:nlev,:) = mmf_input%zmid(:,:,2:nlev,:) - 0.5*abs( &
            mmf_input%zmid(:,:,2:nlev,:) - mmf_input%zmid(:,:,1:nlev-1,:) &
         )
      end if

      ! read thermodynamics
      call nc_get_data(ncid,'T',v4d=mmf_input%t)
      call nc_get_data(ncid,'Q',v4d=mmf_input%sh)
      call nc_get_data(ncid,'RELHUM',v4d=mmf_input%rh)

      ! make sure inputs are from surface to TOA
      if (mmf_input%lev(1) < mmf_input%lev(2)) then
         if (verbose) print *, '   flipping GCM inputs'
         mmf_input%lev = mmf_input%lev(nlev:1:-1)
         mmf_input%pmid = mmf_input%pmid(:,:,nlev:1:-1,:)
         mmf_input%pint = mmf_input%pint(:,:,nlev:1:-1,:)
         mmf_input%zmid = mmf_input%zmid(:,:,nlev:1:-1,:)
         mmf_input%zint = mmf_input%zint(:,:,nlev:1:-1,:)
         mmf_input%t    = mmf_input%t(:,:,nlev:1:-1,:)
         mmf_input%sh   = mmf_input%sh(:,:,nlev:1:-1,:)
         mmf_input%rh   = mmf_input%rh(:,:,nlev:1:-1,:)
         gcm_flipped = .true.
      end if

      ! get CRM inputs
      ! NOTE: the CRM grid has a different number of levels than the GCM grid,
      !       but according to Kharoutinov and Randall (2001), the CRM grid
      !       corresponds to the BOTTOM NZ levels in the GCM grid. We already
      !       flipped the GCM height grid above, so we just need to fill in
      !       the first nz levels with what we read from the CRM, and let those
      !       levels above nz remain zero (we initialized to zero). Of course,
      !       if we had taken care of this before reading in the data such that
      !       nlev == nz, this section of code should still work as expected.
      ! TODO: instead of just filling nz levels here, maybe we should cut off
      !       the CRM grid above, so that everything just has nz levels
      !       consistent with the CRM height grid?
      allocate(v6d(nlon,nlat,nx,ny,nz,ntime))
      call nc_get_data(ncid, 'CRM_QC', v6d=v6d)
      call collapse_xy(v6d, mmf_input%qc(:,:,:,:nz,:))
      call nc_get_data(ncid, 'CRM_QI', v6d=v6d)
      call collapse_xy(v6d, mmf_input%qi(:,:,:,:nz,:))
      call nc_get_data(ncid, 'CRM_QPC', v6d=v6d)
      call collapse_xy(v6d, mmf_input%qpc(:,:,:,:nz,:))
      call nc_get_data(ncid,'CRM_QPI',v6d=v6d)
      call collapse_xy(v6d, mmf_input%qpi(:,:,:,:nz,:))

      ! BRH edit Nov 2015: make sure CRM grid goes from SFC to TOA as well!
      if (mmf_input%crm_z(1) < mmf_input%crm_z(2)) then
         if (verbose) print *, '   flipping CRM inputs'
         mmf_input%crm_z(:) = mmf_input%crm_z(::-1)
         mmf_input%qc(:,:,:,:,:) = mmf_input%qc(:,:,:,::-1,:)
         mmf_input%qi(:,:,:,:,:) = mmf_input%qi(:,:,:,::-1,:)
         mmf_input%qpc(:,:,:,:,:) = mmf_input%qpc(:,:,:,::-1,:)
         mmf_input%qpi(:,:,:,:,:) = mmf_input%qpi(:,:,:,::-1,:)
         crm_flipped = .true.
      end if

      ! graupel may or may not be present; if it is not, we need to calculate
      ! it based on a temperature threshold, stolen from Roj's original driver
      if (nc_var_exists(ncid,'CRM_QG')) then
         call nc_get_data(ncid,'CRM_QG',v6d=v6d)
         call collapse_xy(v6d,mmf_input%qpg(:,:,:,:nz,:))
         mmf_input%qps = mmf_input%qpi
         if (crm_flipped) then
            mmf_input%qpg(:,:,:,:,:) = mmf_input%qpg(:,:,:,::-1,:)
         end if
      else
         allocate(wg_ratio(nlon,nlat,nx*ny,nz,ntime))
         do ix = 1,nx*ny
            wg_ratio(:,:,ix,:nz,:) = (mmf_input%t(:,:,:nz,:) - 223.16) / (283.16-223.16) ! BRH fixed
         end do
         where (wg_ratio<0)
            wg_ratio = 0
         end where
         where (wg_ratio>1)
            wg_ratio = 1
         end where
         do ix = 1,nx*ny
            mmf_input%qpg(:,:,ix,:nz,:) = wg_ratio(:,:,ix,:nz,:)*mmf_input%qpi(:,:,ix,:nz,:)
            mmf_input%qps(:,:,ix,:nz,:) = (1-wg_ratio(:,:,ix,:nz,:))*mmf_input%qpi(:,:,ix,:nz,:)
         end do
         deallocate(wg_ratio)
      end if
      deallocate(v6d)

      ! check condensate mixing ratio units and convert to kg/kg if necessary
      call nc_check(nf90_inq_varid(ncid,'CRM_QC',varid))
      call nc_check(nf90_get_att(ncid,varid,'units',attval))
      if (trim(attval) == 'kg/kg') then
         if (verbose) print *,'  Convert condensate mixing ratios in kg/kg -- leave in kg/kg'
      else if (trim(attval) == 'g/kg') then
         if (verbose) print *,'  Condensate mixing ratios in g/kg -- convert to kg/kg'
         mmf_input%qc = 1e-3*mmf_input%qc
         mmf_input%qi = 1e-3*mmf_input%qi
         mmf_input%qpc = 1e-3*mmf_input%qpc
         mmf_input%qpi = 1e-3*mmf_input%qpi
         mmf_input%qpg = 1e-3*mmf_input%qpg
         mmf_input%qps = 1e-3*mmf_input%qps
      else
         print *,'Condensate mixing ratio units not recognized.'
         stop
      end if

      ! check heights and make sure consistent
      if (mmf_input%crm_z(1) .ne. mmf_input%lev(1)) then
         print *, 'GCM and CRM heights do not match.'
         stop
      end if
      if ((mmf_input%crm_z(1) < mmf_input%crm_z(2)) &
            .or. (mmf_input%lev(1) < mmf_input%lev(2)) &
            .or. (mmf_input%pmid(1,1,1,1) < mmf_input%pmid(1,1,2,1)) &
            .or. (mmf_input%zmid(1,1,1,1) > mmf_input%zmid(1,1,2,1))) then
         print *, 'Inconsistent heights.'
         stop
      end if

      ! close file
      call nc_check(nf90_close(ncid))
   end subroutine read_mmf_input


   ! ---------------------------------------------------------------------------
   ! subroutine to combine crm_x and crm_y dimensions
   ! ---------------------------------------------------------------------------
   subroutine collapse_xy(v6d, v5d)
      real,intent(in) :: v6d(:,:,:,:,:,:)
      real,intent(inout) :: v5d(:,:,:,:,:)
      integer :: nx, ny, ix, iy, m

      nx = size(v6d, 3)
      ny = size(v6d, 4)

      m = 1
      do iy = 1, ny
         do ix = 1, nx
            v5d(:,:,m,:,:) = v6d(:,:,ix,iy,:,:)
            m = m + 1
         end do
      end do
   end subroutine


   ! ---------------------------------------------------------------------------
   ! subroutine to populate COSP gridbox from MMF fields
   ! ---------------------------------------------------------------------------
   subroutine populate_gridbox(mmf_input, gbx, lon, lat, time, optics_flag)
      type(mmf_input_type), intent(in) :: mmf_input
      type(cosp_gridbox), intent(inout) :: gbx
      integer,  intent(in) :: lon, lat, time
      integer, optional, intent(in) :: optics_flag
         ! how to calculate cloud optics
         ! 1 -> use CAM subroutines
         ! 2 -> use RRTMG subroutines
         ! 3 -> use ISCCP approximations

      ! local variables
      integer :: m, nlev, ncol
      integer :: optics_flag_local = 0

      ! get options
      if (present(optics_flag)) optics_flag_local = optics_flag

      ! get dimension sizes
      if (mmf_input%azg_input) then
         nlev = mmf_input%nzm_used(lon,lat,time)
      else
         nlev = mmf_input%nlev
      end if

      ! large-scale fields from MMF GCM
      gbx%longitude(:) = mmf_input%lon(lon)
      gbx%latitude(:) = mmf_input%lat(lat)
      gbx%psfc(:) = mmf_input%psfc(lon,lat,time)
      gbx%skt(:) = mmf_input%tsfc(lon,lat,time)
      gbx%sunlit(:) = mmf_input%sunlit(lon,lat,time)
      do m = 1,gbx%npoints
         gbx%p(m,:) = mmf_input%pmid(lon,lat,:nlev,time)
         gbx%ph(m,:) = mmf_input%pint(lon,lat,:nlev,time)
         gbx%zlev(m,:) = mmf_input%zmid(lon,lat,:nlev,time)
         gbx%zlev_half(m,:) = mmf_input%zint(lon,lat,:nlev,time)
         gbx%t(m,:) = mmf_input%t(lon,lat,:nlev,time)
         gbx%sh(m,:) = mmf_input%sh(lon,lat,:nlev,time)
         gbx%q(m,:) = mmf_input%rh(lon,lat,:nlev,time)
      end do ! loop over points
      where (gbx%zlev_half<0)
         gbx%zlev_half = 0
      end where

      ! hydrometeor fields from MMF CRM
      if (gbx%ncolumns == 1) then ! CRM mode
         gbx%mr_hydro(:,:nlev,I_LSCLIQ) = mmf_input%qc(lon,lat,:,:nlev,time)
         gbx%mr_hydro(:,:nlev,I_LSCICE) = mmf_input%qi(lon,lat,:,:nlev,time)
         gbx%mr_hydro(:,:nlev,I_CVCLIQ) = 0.0
         gbx%mr_hydro(:,:nlev,I_CVCICE) = 0.0
         gbx%mr_hydro(:,:nlev,I_LSRAIN) = mmf_input%qpc(lon,lat,:,:nlev,time)
         gbx%mr_hydro(:,:nlev,I_LSSNOW) = mmf_input%qps(lon,lat,:,:nlev,time)
         gbx%mr_hydro(:,:nlev,I_LSGRPL) = mmf_input%qpg(lon,lat,:,:nlev,time)
         gbx%mr_hydro(:,:nlev,I_CVRAIN) = 0.0
         gbx%mr_hydro(:,:nlev,I_CVSNOW) = 0.0
         gbx%np(:,:nlev,I_LSCLIQ) = mmf_input%nc(lon,lat,:,:nlev,time)
         gbx%np(:,:nlev,I_LSCICE) = mmf_input%ni(lon,lat,:,:nlev,time)
         gbx%np(:,:nlev,I_LSRAIN) = mmf_input%nr(lon,lat,:,:nlev,time)
         gbx%np(:,:nlev,I_LSSNOW) = mmf_input%ns(lon,lat,:,:nlev,time)
         gbx%np(:,:nlev,I_LSGRPL) = mmf_input%ng(lon,lat,:,:nlev,time)
      else ! GCM mode; pass gridbox mean values
         ncol = size(mmf_input%qc,3)
         do m = 1,gbx%npoints
            gbx%mr_hydro(m,:nlev,I_LSCLIQ) = sum(mmf_input%qc(lon,lat,:,:nlev,time),1)/ncol
            gbx%mr_hydro(m,:nlev,I_LSCICE) = sum(mmf_input%qi(lon,lat,:,:nlev,time),1)/ncol
            gbx%mr_hydro(m,:nlev,I_CVCLIQ) = 0.0
            gbx%mr_hydro(m,:nlev,I_CVCICE) = 0.0
            gbx%mr_hydro(m,:nlev,I_LSRAIN) = sum(mmf_input%qpc(lon,lat,:,:nlev,time),1)/ncol
            gbx%mr_hydro(m,:nlev,I_LSSNOW) = sum(mmf_input%qps(lon,lat,:,:nlev,time),1)/ncol
            gbx%mr_hydro(m,:nlev,I_LSGRPL) = sum(mmf_input%qpg(lon,lat,:,:nlev,time),1)/ncol
            gbx%mr_hydro(m,:nlev,I_CVRAIN) = 0.0
            gbx%mr_hydro(m,:nlev,I_CVSNOW) = 0.0
            gbx%np(m,:nlev,I_LSCLIQ) = sum(mmf_input%nc(lon,lat,:,:nlev,time),1)/ncol
            gbx%np(m,:nlev,I_LSCICE) = sum(mmf_input%ni(lon,lat,:,:nlev,time),1)/ncol
            gbx%np(m,:nlev,I_LSRAIN) = sum(mmf_input%nr(lon,lat,:,:nlev,time),1)/ncol
            gbx%np(m,:nlev,I_LSSNOW) = sum(mmf_input%ns(lon,lat,:,:nlev,time),1)/ncol
            gbx%np(m,:nlev,I_LSGRPL) = sum(mmf_input%ng(lon,lat,:,:nlev,time))/ncol
         end do ! loop over points
      end if ! pass full CRM inputs or averages for GCM mode

      ! cloud optical properties
      if (optics_flag_local == 0) then ! use CAM cloud optics
         do m = 1,gbx%npoints
            call cam_optics( &
               gbx%psfc(m), &
               gbx%p(m,:), &
               gbx%ph(m,:), &
               gbx%t(m,:), &
               gbx%mr_hydro(m,:,I_LSCLIQ), &
               gbx%mr_hydro(m,:,I_LSCICE), &
               gbx%reff(m,:,I_LSCLIQ), &
               gbx%reff(m,:,I_LSCICE), &
               gbx%dtau_s(m,:), &
               gbx%dem_s(m,:), &
               gbx%tca(m,:) &
            )
         end do ! loop over columns
      else if (optics_flag_local == 1) then ! use RRTMG cloud optics
         print *, 'RRTMG cloud optics package not supported.'
         stop
      else if (optics_flag_local == 2) then ! use ISCCP cloud optics approximations
         print *, 'ISCCP cloud optics package not supported.'
         stop
      end if ! cloud optics calculation type
      gbx%dtau_c(:,:nlev) = 0.0  ! no convective cloud designation in MMF
      gbx%dem_c(:,:nlev) = 0.0   ! no convective cloud designation in MMF
      gbx%cca(:,:nlev) = 0.0     ! no convective cloud designation in MMF
      gbx%reff = 0.0             ! do not use effective radii
   end subroutine populate_gridbox


   ! ---------------------------------------------------------------------------
   ! subroutine to calculate cloud optical depth
   ! (code copied from CAM routine radcswmx in radsw.F90)
   ! abarl: A coefficient for extinction optical depth
   ! abari: A coefficient for extinction optical depth
   ! bbarl: B coefficient for extinction optical depth
   ! bbari: B coefficient for extinction optical depth
   ! ---------------------------------------------------------------------------
   subroutine cldtau(cliqwp, cicewp, rel, rei, cld, ctau)
      real, intent(in) :: cliqwp(:), cicewp(:)
      real, intent(in) :: rel(:), rei(:)
      real, intent(in) :: cld(:)
      real, intent(inout) :: ctau(:)
      real, parameter :: abarli = 2.817e-02
      real, parameter :: bbarli = 1.305
      real, parameter :: abarii = 3.448e-03
      real, parameter :: bbarii = 2.431
      real, parameter :: cldmin = 1e-30
      real, parameter :: cldeps = 1e-30
      real :: tmp1l, tmp1i
      integer :: i, k, pver
      real, allocatable :: tauxcl(:), tauxci(:)

      pver = size(cliqwp, 1)
      allocate(tauxcl(pver), tauxci(pver))
      do k = 1, pver
         if (cld(k) >= cldmin .and. cld(k) >= cldeps) then
            tmp1l = abarli + bbarli/rel(k)
            tmp1i = abarii + bbarii/rei(k)
            tauxcl(k) = cliqwp(k)*tmp1l
            tauxci(k) = cicewp(k)*tmp1i
         else
            tauxcl(k) = 0.0
            tauxci(k) = 0.0
         end if
         ctau(k) = tauxcl(k) + tauxci(k)
      end do
      deallocate(tauxcl, tauxci)
   end subroutine cldtau


   ! ---------------------------------------------------------------------------
   ! subroutine to copy COSP point outputs to global structure
   ! ---------------------------------------------------------------------------
   subroutine save_output(cfg, &
                          gridbox1, subgrid1, sgradar1, sglidar1, &
                          isccp1, misr1, modis1, stradar1, stlidar1, &
                          gridbox2, subgrid2, sgradar2, sglidar2,  &
                          isccp2, misr2, modis2, stradar2, stlidar2, loc)
      type(cosp_config), intent(in) :: cfg
      type(cosp_gridbox), intent(in) :: gridbox1
      type(cosp_subgrid), intent(in) :: subgrid1
      type(cosp_sgradar), intent(in) :: sgradar1
      type(cosp_sglidar), intent(in) :: sglidar1
      type(cosp_isccp), intent(in) :: isccp1
      type(cosp_misr), intent(in) :: misr1
      type(cosp_modis), intent(in) :: modis1
      type(cosp_radarstats), intent(in) :: stradar1
      type(cosp_lidarstats), intent(in) :: stlidar1
      type(cosp_gridbox), intent(inout) :: gridbox2
      type(cosp_subgrid), intent(inout) :: subgrid2
      type(cosp_sgradar), intent(inout) :: sgradar2
      type(cosp_sglidar), intent(inout) :: sglidar2
      type(cosp_isccp), intent(inout) :: isccp2
      type(cosp_misr), intent(inout) :: misr2
      type(cosp_modis), intent(inout) :: modis2
      type(cosp_radarstats), intent(inout) :: stradar2
      type(cosp_lidarstats), intent(inout) :: stlidar2
      integer, intent(in), optional :: loc
      integer :: npoints, ncolumns

      npoints = sgradar1%npoints
      ncolumns = sgradar1%ncolumns
      
      ! gridbox coordinates
      gridbox2%latitude(loc) = gridbox1%latitude(1)
      gridbox2%longitude(loc) = gridbox1%longitude(1)

      ! subgrid diagnostics
      if (cfg%Lisccp_sim) then
         if (ncolumns == 1) then
            subgrid2%frac_out(loc,:,:) = subgrid1%frac_out(:,1,:)
         else
            subgrid2%prec_frac(loc,:,:) = subgrid1%prec_frac(1,:,:)
            subgrid2%frac_out(loc,:,:) = subgrid1%frac_out(1,:,:)
         end if
      end if

      ! isccp diagnostics
      if (cfg%Lisccp_sim) then
         if (ncolumns == 1) then
            isccp2%boxtau(loc,:) = isccp1%boxtau(:,1)
            isccp2%boxptop(loc,:) = isccp1%boxptop(:,1)
         else
            isccp2%boxtau(loc,:) = isccp1%boxtau(1,:)
            isccp2%boxptop(loc,:) = isccp1%boxptop(1,:)
         end if

         isccp2%fq_isccp(loc,:,:) = 1.0*sum( &
            isccp1%fq_isccp, 1 &
         )/size(isccp1%fq_isccp, 1)

         ! only average over cloudy part of gridbox
         isccp2%meanptop(loc) = 1.0*sum( &
            isccp1%meanptop, &
            mask = (isccp1%meanptop > 0) &
         )/max(1, count(isccp1%meanptop > 0))

         isccp2%meantaucld(loc) = 1.0*sum( &
            isccp1%meantaucld, &
            mask = (isccp1%meanptop > 0) &
         )/max(1, count(isccp1%meanptop > 0))

         isccp2%meanalbedocld(loc) = 1.0*sum( &
            isccp1%meanalbedocld,  &
            mask = (isccp1%meanptop > 0) &
         )/max(1, count(isccp1%meanptop > 0))

         isccp2%meantb(loc) = 1.0*sum( &
            isccp1%meantb,  &
            mask = (isccp1%meanptop > 0) &
         )/max(1, count(isccp1%meanptop > 0))

         isccp2%meantbclr(loc) = 1.0*sum( &
            isccp1%meantbclr,  &
            mask = (isccp1%meanptop > 0) &
         )/max(1, count(isccp1%meanptop > 0))

         ! total cloud fraction...average over all points
         isccp2%totalcldarea(loc) = sum( &
            isccp1%totalcldarea &
         )/size(isccp1%totalcldarea)

      end if

      ! modis diagnostics
      if (cfg%Lmodis_sim) then
         modis2%optical_thickness_vs_cloud_top_pressure(loc,:,:) = 1.0*sum( &
            modis1%optical_thickness_vs_cloud_top_pressure, 1 &
         )/size(modis1%optical_thickness_vs_cloud_top_pressure, 1)

         ! only average over cloudy part of gridbox
         modis2%cloud_top_pressure_total_mean(loc) = 1.0*sum( &
            modis1%cloud_top_pressure_total_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%ice_water_path_mean(loc) = 1.0*sum( &
            modis1%ice_water_path_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%liquid_water_path_mean(loc) = 1.0*sum( &
            modis1%liquid_water_path_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%cloud_particle_size_ice_mean(loc) = 1.0*sum( &
            modis1%cloud_particle_size_ice_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%cloud_particle_size_water_mean(loc) = 1.0*sum( &
            modis1%cloud_particle_size_water_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_ice_logmean(loc) = 1.0*sum( &
            modis1%optical_thickness_ice_logmean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_water_logmean(loc) = 1.0*sum( &
            modis1%optical_thickness_water_logmean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_total_logmean(loc) = 1.0*sum( &
            modis1%optical_thickness_total_logmean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_ice_mean(loc) = 1.0*sum( &
            modis1%optical_thickness_ice_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_water_mean(loc) = 1.0*sum( &
            modis1%optical_thickness_water_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         modis2%optical_thickness_total_mean(loc) = 1.0*sum( &
            modis1%optical_thickness_total_mean, &
            mask=(modis1%cloud_top_pressure_total_mean > 0) &
         )/max(1, count(modis1%cloud_top_pressure_total_mean > 0))

         ! modis cloud fractions...average over all points
         modis2%cloud_fraction_low_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_low_mean &
         )/size(modis1%cloud_fraction_low_mean)

         modis2%cloud_fraction_mid_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_mid_mean &
         )/size(modis1%cloud_fraction_mid_mean)

         modis2%cloud_fraction_high_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_high_mean &
         )/size(modis1%cloud_fraction_high_mean)

         modis2%cloud_fraction_ice_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_ice_mean &
         )/size(modis1%cloud_fraction_ice_mean)

         modis2%cloud_fraction_water_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_water_mean &
         )/size(modis1%cloud_fraction_water_mean)

         modis2%cloud_fraction_total_mean(loc) = 1.0*sum( &
            modis1%cloud_fraction_total_mean &
         )/size(modis1%cloud_fraction_total_mean)
      end if ! do modis

      ! misr diagnostics
      if (cfg%Lmisr_sim) then
         misr2%fq_misr(loc,:,:) = 1.0*sum( &
            misr1%fq_misr, 1 &
         )/size(misr1%fq_misr, 1)
      end if ! do misr

      ! radar diagnostics
      if (cfg%Lradar_sim) then
         if (ncolumns == 1) then
            sgradar2%ze_tot(loc,:,:) = sgradar1%ze_tot(:,1,:)
         else
            sgradar2%ze_tot(loc,:,:) = sgradar1%ze_tot(1,:,:)
         end if
         sgradar2%att_gas(loc,:) = sgradar1%att_gas(1,:)

         stradar2%cfad_ze(loc,:,:) = sum( &
            stradar1%cfad_ze, 1 &
         )/size(stradar1%cfad_ze, 1)

         stradar2%lidar_only_freq_cloud(loc,:) = sum( &
            stradar1%lidar_only_freq_cloud, 1 &
         )/size(stradar1%lidar_only_freq_cloud, 1)

         stradar2%radar_lidar_tcc(loc) = sum( &
            stradar1%radar_lidar_tcc, 1 &
         )/size(stradar1%radar_lidar_tcc, 1)
      end if ! do radar

      ! lidar diagnostics
      if (cfg%Llidar_sim) then
         if (ncolumns == 1) then
            sglidar2%beta_tot(loc,:,:) = sglidar1%beta_tot(:,1,:)
            sglidar2%tau_tot(loc,:,:) = sglidar1%tau_tot(:,1,:)
            sglidar2%refl(loc,:,:) = sglidar1%refl(:,1,:)
         else
            sglidar2%beta_tot(loc,:,:) = sglidar1%beta_tot(1,:,:)
            sglidar2%tau_tot(loc,:,:) = sglidar1%tau_tot(1,:,:)
            sglidar2%refl(loc,:,:) = sglidar1%refl(1,:,:)
         end if
         sglidar2%beta_mol(loc,:) = sglidar1%beta_mol(1,:)
         stlidar2%srbval = stlidar1%srbval

         stlidar2%cfad_sr(loc,:,:) = sum( &
            stlidar1%cfad_sr, 1 &
         )/size(stlidar1%cfad_sr, 1)

         stlidar2%parasolrefl(loc,:) = sum( &
            stlidar1%parasolrefl, 1 &
         )/size(stlidar1%parasolrefl, 1)

         ! Roj sums all points greater than or equal to zero, but
         ! counts all points whether undefined or not, effectively
         ! treating undefined as "no cloud".
         stlidar2%cldlayer(loc,:) = sum( &
            stlidar1%cldlayer, 1, &
            mask=(stlidar1%cldlayer >= 0) &
         )/size(stlidar1%cldlayer, 1)

         stlidar2%lidarcld(loc,:) = sum( &
            stlidar1%lidarcld, 1, &
            mask=(stlidar1%lidarcld >= 0) &
         )/size(stlidar1%lidarcld, 1)
      end if ! do lidar
   end subroutine save_output


   subroutine cam_optics(psfc, p, ph, t, qc, qi, rel, rei, dtau, dem, tca)
      real, intent(in) :: psfc
      real, intent(in) :: p(:), ph(:), t(:)
      real, intent(in) :: qc(:), qi(:)
      real, intent(inout) :: rel(:), rei(:)
      real, intent(inout) :: tca(:), dtau(:), dem(:)

      real, allocatable :: clwp(:) ! cloud liquid water path
      real, allocatable :: ciwp(:) ! cloud ice water path
      real, allocatable :: ctwp(:) ! cloud total water path
      real, allocatable :: fice(:) ! cloud ice fraction
      real(r8), allocatable :: rel_r8(:) ! liquid water effective radius
      real(r8), allocatable :: rei_r8(:) ! ice effective radius
      real(r8), allocatable :: dem_r8(:) ! later emissivity
      integer :: nlevels

      ! allocate memory for temporary arrays
      nlevels = size(qc, 1)
      allocate(clwp(nlevels), ciwp(nlevels), ctwp(nlevels), fice(nlevels), &
               rel_r8(nlevels), rei_r8(nlevels), dem_r8(nlevels))

      ! calculate cloud water paths
      where (qc(:) + qi(:) > 1e-30)
         clwp = qc(:) * 2.0 * (ph(:) - p(:)) / 9.81 * 1e3
         ciwp = qi(:) * 2.0 * (ph(:) - p(:)) / 9.81 * 1e3
      elsewhere
         clwp = 0
         ciwp = 0
      end where
      ctwp = clwp + ciwp

      ! calculate ice ratio
      where (ctwp > 1e-30)
         tca(:) = 1
         fice = ciwp/ctwp
      elsewhere
         tca(:) = 0
         fice = 0
      end where

      ! calculate ice ratio
      where (ctwp > 1e-30)
         tca(:) = 1
         fice = ciwp/ctwp
      elsewhere
         tca(:) = 0
         fice = 0
      end where

      ! calculate effective radii
      call cldefr(1, 1, 1, &
                  (/0.0_r8/), t(:)*1.0_r8, &
                  rel_r8(:), rei_r8(:), (/psfc*1.0_r8/), p(:)*1.0_r8, &
                  (/0.0_r8/), (/0.0_r8/), (/0.0_r8/), nlevels)

      ! calculate emissivity
      call cldems(1, 1, 1, &
                  ctwp(:)*1.0_r8, fice(:)*1.0_r8, &
                  rei_r8(:), dem_r8(:), nlevels)

      ! copy data from explicit kind variables
      rel(:) = rel_r8(:)
      rei(:) = rei_r8(:)
      dem(:) = dem_r8(:)

      ! calculate cloud optical depth
      call cldtau(clwp(:), ciwp(:), &
                  rel(:), rei(:), &
                  tca(:), dtau(:))

      ! convert from micron (COSP assumes effective radii in meters)
      rel(:) = rel(:)*1e-6
      rei(:) = rei(:)*1e-6

      ! free memory
      deallocate(clwp, ciwp, ctwp, fice, rel_r8, rei_r8, dem_r8)

   end subroutine cam_optics


   !
   ! DEPRECATED
   !
   subroutine write_mmf_output(filename,mmf_input)
      character(len=*),intent(in) :: filename 
      type(mmf_input_type),intent(in) :: mmf_input
      integer :: ncid,varid,dimid
      integer :: nlon,nlat,ncol,nlev,ntime
      integer :: time_id,lat_id,lon_id,col_id,lev_id

      ! get dimension sizes
      nlon = size(mmf_input%qc,1)
      nlat = size(mmf_input%qc,2)
      ncol = size(mmf_input%qc,3)
      nlev = size(mmf_input%qc,4)
      ntime = size(mmf_input%qc,5)

      ! open file
      call nc_check(nf90_create(trim(filename),NF90_SHARE,ncid))

      ! define dimensions
      call nc_check(nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_id))
      call nc_check(nf90_def_dim(ncid,'lat',nlat,lat_id))
      call nc_check(nf90_def_dim(ncid,'lon',nlon,lon_id))
      call nc_check(nf90_def_dim(ncid,'col',ncol,col_id))
      call nc_check(nf90_def_dim(ncid,'lev',nlev,lev_id))

      ! define variables
      call nc_check(nf90_def_var(ncid,'time',NF90_FLOAT,(/time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','time'))
      call nc_check(nf90_put_att(ncid,varid,'units','days since 1998-01-01 00:00:00'))

      call nc_check(nf90_def_var(ncid,'lat',NF90_FLOAT,(/lat_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','latitude'))
      call nc_check(nf90_put_att(ncid,varid,'units','degrees_north'))

      call nc_check(nf90_def_var(ncid,'lon',NF90_FLOAT,(/lon_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','longitude'))
      call nc_check(nf90_put_att(ncid,varid,'units','degrees_east'))

      call nc_check(nf90_def_var(ncid,'psfc',NF90_FLOAT,(/lon_id,lat_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','surface pressure'))
      call nc_check(nf90_put_att(ncid,varid,'units','Pa'))

      call nc_check(nf90_def_var(ncid,'tsfc',NF90_FLOAT,(/lon_id,lat_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','surface temperature'))
      call nc_check(nf90_put_att(ncid,varid,'units','K'))

      call nc_check(nf90_def_var(ncid,'sunlit',NF90_FLOAT,(/lon_id,lat_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','sunlit flag'))
      call nc_check(nf90_put_att(ncid,varid,'units','1'))

      call nc_check(nf90_def_var(ncid,'pmid',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','pressure at midpoints'))
      call nc_check(nf90_put_att(ncid,varid,'units','Pa'))

      call nc_check(nf90_def_var(ncid,'pint',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','pressure at interfaces'))
      call nc_check(nf90_put_att(ncid,varid,'units','Pa'))

      call nc_check(nf90_def_var(ncid,'zmid',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','height at midpoints'))
      call nc_check(nf90_put_att(ncid,varid,'units','m'))

      call nc_check(nf90_def_var(ncid,'zint',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','height at interfaces'))
      call nc_check(nf90_put_att(ncid,varid,'units','m'))

      call nc_check(nf90_def_var(ncid,'t',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','temperature'))
      call nc_check(nf90_put_att(ncid,varid,'units','K'))

      call nc_check(nf90_def_var(ncid,'sh',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','water vapor specific humidity'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      call nc_check(nf90_def_var(ncid,'rh',NF90_FLOAT,(/lon_id,lat_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','relative humidity'))
      call nc_check(nf90_put_att(ncid,varid,'units','%'))

      call nc_check(nf90_def_var(ncid,'qc',NF90_FLOAT,(/lon_id,lat_id,col_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','cloud liquid mixing ratio'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      call nc_check(nf90_def_var(ncid,'qi',NF90_FLOAT,(/lon_id,lat_id,col_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','cloud ice mixing ratio'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      call nc_check(nf90_def_var(ncid,'qpc',NF90_FLOAT,(/lon_id,lat_id,col_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','precipitating rain mixing ratio'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      call nc_check(nf90_def_var(ncid,'qps',NF90_FLOAT,(/lon_id,lat_id,col_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','precipitating snow mixing ratio'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      call nc_check(nf90_def_var(ncid,'qpg',NF90_FLOAT,(/lon_id,lat_id,col_id,lev_id,time_id/),varid))
      call nc_check(nf90_put_att(ncid,varid,'long_name','precipitating graupel mixing ratio'))
      call nc_check(nf90_put_att(ncid,varid,'units','kg/kg'))

      ! end define mode
      call nc_check(nf90_enddef(ncid))

      !
      ! write coordinate variables
      !
      call nc_check(nf90_inq_varid(ncid,'time',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%time))

      call nc_check(nf90_inq_varid(ncid,'lat',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%lat))

      call nc_check(nf90_inq_varid(ncid,'lon',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%lon))

      !
      ! write variables
      !
      call nc_check(nf90_inq_varid(ncid,'psfc',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%psfc))

      call nc_check(nf90_inq_varid(ncid,'tsfc',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%tsfc))

      call nc_check(nf90_inq_varid(ncid,'sunlit',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%sunlit))

      call nc_check(nf90_inq_varid(ncid,'pmid',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%pmid))

      call nc_check(nf90_inq_varid(ncid,'pint',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%pint))

      call nc_check(nf90_inq_varid(ncid,'zmid',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%zmid))

      call nc_check(nf90_inq_varid(ncid,'zint',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%zint))

      call nc_check(nf90_inq_varid(ncid,'t',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%t))

      call nc_check(nf90_inq_varid(ncid,'sh',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%sh))

      call nc_check(nf90_inq_varid(ncid,'rh',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%rh))

      call nc_check(nf90_inq_varid(ncid,'qc',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%qc))

      call nc_check(nf90_inq_varid(ncid,'qi',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%qi))

      call nc_check(nf90_inq_varid(ncid,'qpc',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%qpc))

      call nc_check(nf90_inq_varid(ncid,'qps',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%qps))

      call nc_check(nf90_inq_varid(ncid,'qpg',varid))
      call nc_check(nf90_put_var(ncid,varid,mmf_input%qpg))

      ! close files
      call nc_check(nf90_close(ncid))

   end subroutine write_mmf_output

end module cosp_mmf_utils
