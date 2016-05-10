      module netcdf_utils
         use netcdf

         private
         public nc_check
         public nc_var_exists
         public nc_get_dim
         public nc_get_data


      contains
         ! ---------------------------------------------------------------------------
         ! subroutine to handle NetCDF errors
         ! ---------------------------------------------------------------------------
         subroutine nc_check(stat)
            integer,intent(in) :: stat
            if (stat .ne. nf90_noerr) then
               print *, trim(nf90_strerror(stat))
               stop 'stopped'
            end if
         end subroutine nc_check

         ! ---------------------------------------------------------------------------
         ! function to check if variable exists in NetCDF file
         ! ---------------------------------------------------------------------------
         logical function nc_var_exists(ncid,vname)
            integer,intent(in) :: ncid
            character(len=*),intent(in) :: vname
            integer :: varid,nvars
            character(len=32) :: ncname

            nc_var_exists = .false.
            call nc_check(nf90_inquire(ncid,nVariables=nvars))
            do varid = 1,nvars
               call nc_check(nf90_inquire_variable(ncid,varid,name=ncname))
               if (trim(ncname).eq.trim(vname)) then
                  nc_var_exists = .true.
               end if
            end do
            return
         end function nc_var_exists

         ! ---------------------------------------------------------------------------
         ! subroutine to get dimension lengths from NetCDF file by name
         ! ---------------------------------------------------------------------------
         subroutine nc_get_dim(ncid,dname,dlen)
            integer,intent(in) :: ncid
            character(len=*),intent(in) :: dname
            integer,intent(out) :: dlen
            integer :: dimid, stat

            stat=nf90_inq_dimid(ncid,dname,dimid)
            if (stat .ne. nf90_noerr) then
               print *, '   Can not find dimension name: ',dname
               dlen = -1   
            else
               call nc_check(nf90_inquire_dimension(ncid,dimid,len=dlen))
            end if
         end subroutine nc_get_dim

         ! ---------------------------------------------------------------------------
         ! subroutine to read real or double precision variables by name
         ! ---------------------------------------------------------------------------
     
         subroutine nc_get_data( &
            ncid,vname, &
            v1d,v2d,v3d,v4d,v5d,v6d, &
            d1d,d2d,d3d,d4d,d5d,d6d &
         )
            integer,intent(in) :: ncid
            character(len=*),intent(in) :: vname
            real,intent(inout),optional :: v1d(:)
            real,intent(inout),optional :: v2d(:,:)
            real,intent(inout),optional :: v3d(:,:,:)
            real,intent(inout),optional :: v4d(:,:,:,:)
            real,intent(inout),optional :: v5d(:,:,:,:,:)
            real,intent(inout),optional :: v6d(:,:,:,:,:,:)
            double precision,intent(inout),optional :: d1d(:)
            double precision,intent(inout),optional :: d2d(:,:)
            double precision,intent(inout),optional :: d3d(:,:,:)
            double precision,intent(inout),optional :: d4d(:,:,:,:)
            double precision,intent(inout),optional :: d5d(:,:,:,:,:)
            double precision,intent(inout),optional :: d6d(:,:,:,:,:,:)

            integer :: varid

            call nc_check(nf90_inq_varid(ncid,vname,varid))
            if (present(v1d)) then
               call nc_check(nf90_get_var(ncid,varid,v1d))
            else if (present(v2d)) then
               call nc_check(nf90_get_var(ncid,varid,v2d))
            else if (present(v3d)) then
               call nc_check(nf90_get_var(ncid,varid,v3d))
            else if (present(v4d)) then
               call nc_check(nf90_get_var(ncid,varid,v4d))
            else if (present(v5d)) then
               call nc_check(nf90_get_var(ncid,varid,v5d))
            else if (present(v6d)) then
               call nc_check(nf90_get_var(ncid,varid,v6d))
            end if
            if (present(d1d)) then
               call nc_check(nf90_get_var(ncid,varid,d1d))
            else if (present(d2d)) then
               call nc_check(nf90_get_var(ncid,varid,d2d))
            else if (present(d3d)) then
               call nc_check(nf90_get_var(ncid,varid,d3d))
            else if (present(d4d)) then
               call nc_check(nf90_get_var(ncid,varid,d4d))
            else if (present(d5d)) then
               call nc_check(nf90_get_var(ncid,varid,d5d))
            else if (present(d6d)) then
               call nc_check(nf90_get_var(ncid,varid,d6d))
            end if
         end subroutine nc_get_data

end module
