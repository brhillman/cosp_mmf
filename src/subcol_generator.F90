module subcol_generator

implicit none

public gen_subcol_cld,gen_subcol_precip,gen_subcol_prec_scops

contains
   ! Implementation of Raisanen et al. 2004 subcolumn generator
   subroutine gen_subcol_cld(pmid,cf,qc_avg,qi_avg,qc_col,qi_col)
      real,intent(in) :: pmid(:) ! used to seed random number generator
      real,intent(in) :: cf(:) ! cloud fraction
      real,intent(in) :: qc_avg(:) ! average (in-cloud) liquid water mixing ratio
      real,intent(in) :: qi_avg(:) ! average (in-cloud) ice water mixing ratio
      real,intent(inout) :: qc_col(:,:) ! column cloud liquid water mixing ratio
      real,intent(inout) :: qi_col(:,:) ! column cloud ice water mixing ratio
      
      ! local variables
      real,allocatable :: rnd_vec(:)
      real,allocatable :: rn(:,:)
      real,allocatable :: x(:,:)
      integer :: ncol,nlev,icol,ilev
      integer :: seed(1)

      ! dimension sizes
      ncol = size(qc_col,1)
      nlev = size(qc_col,2)

      ! allocate arrays
      allocate(rnd_vec(ncol))
      allocate(rn(ncol,nlev))
      allocate(x(ncol,nlev))

      ! generate random numbers
      seed = (pmid(1) - int(pmid(1))) * 1000000000
      call random_seed(put=seed)
      call random_number(rn)

      ! calculate x(i,j)
      ! NOTE: inputs must be from TOA to surface
      ! NOTE: this implements max-rand overlap...
      !       need to modify to implement generalized overlap
      do icol = 1,ncol
         x(icol,1) = rn(icol,1)
         do ilev = 2,nlev
            if (x(icol,ilev-1) > (1.0 - cf(ilev-1))) then
               ! maximal overlap; use same random number as previous level
               x(icol,ilev) = x(icol,ilev-1)
            else
               ! random overlap; use new random number
               x(icol,ilev) = rn(icol,ilev)*(1.0 - cf(ilev-1))
            end if
         end do
      end do

      ! cloud properties
      do icol = 1,ncol
         do ilev = 1,nlev
            if (x(icol,ilev) > (1.0 - cf(ilev))) then
               ! cloudy
               qc_col(icol,ilev) = qc_avg(ilev)
               qi_col(icol,ilev) = qi_avg(ilev)
            else
               qc_col(icol,ilev) = 0.0
               qi_col(icol,ilev) = 0.0
            end if
         end do
      end do

      ! deallocate arrays
      deallocate(rnd_vec)
      deallocate(rn)
      deallocate(x)

   end subroutine


!  subroutine gen_subcol_scops(pmid,cf,qc_avg,qi_avg,qc_col,qi_col)
!     real,intent(in) :: pmid(:) ! used to seed random number generator
!     real,intent(in) :: cf(:) ! cloud fraction
!     real,intent(in) :: qc_avg(:) ! average (in-cloud) liquid water mixing ratio
!     real,intent(in) :: qi_avg(:) ! average (in-cloud) ice water mixing ratio
!     real,intent(inout) :: qc_col(:,:) ! column cloud liquid water mixing ratio
!     real,intent(inout) :: qi_col(:,:) ! column cloud ice water mixing ratio
!  end subroutine gen_subcol_scops
 

   subroutine gen_subcol_precip(precip_frac,qpc_avg,qps_avg,qpg_avg,qc_col,qpc_col,qps_col,qpg_col,precip_mode)
      real,intent(in) :: precip_frac(:) ! precipitation fraction
      real,intent(in) :: qpc_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qps_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qpg_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qc_col(:,:) ! column total cloud mixing ratio
      real,intent(inout) :: qpc_col(:,:) ! column precipitation mixing ratio
      real,intent(inout) :: qps_col(:,:) ! column precipitation mixing ratio
      real,intent(inout) :: qpg_col(:,:) ! column precipitation mixing ratio
      integer,intent(in) :: precip_mode

      ! local variables
      integer :: ncol,nlev,icol,ilev
      real,allocatable :: precip_flag(:,:)

      ! dimension sizes
      ncol = size(qc_col,1)
      nlev = size(qc_col,2)

      allocate(precip_flag(ncol,nlev))
      precip_flag = 0

      do icol = 1,ncol
         ! first level (TOA)
         if (precip_frac(1) > 0 .and. qc_col(icol,1) > 0) then
            precip_flag(icol,1) = 1
         else
            precip_flag(icol,1) = 0
         end if

         ! walk down column
         do ilev = 2,nlev
            if (precip_frac(ilev)>0) then
               ! assign precip to cells with cloud in current level or with precip above
               if (qc_col(icol,ilev)>0 .or. precip_flag(icol,ilev-1)==1) then
                  precip_flag(icol,ilev) = 1
               else
                  precip_flag(icol,ilev) = 0
               end if
            else
               precip_flag(icol,ilev) = 0
            end if
         end do
      end do

      ! adjust precipitation
      if (precip_mode == 1) then ! constrain precip by input precip fraction
         call constrain_subcol_precip(precip_frac,precip_flag)
      else if (precip_mode == 2) then ! remove non-cloudy precip
         where(qc_col == 0)
            precip_flag = 0
         endwhere
      end if

      ! assign precip mixing ratios
      do icol = 1,ncol
         where (precip_flag(icol,:) == 1)
            qpc_col(icol,:) = qpc_avg(:)
            qps_col(icol,:) = qps_avg(:)
            qpg_col(icol,:) = qpg_avg(:)
         else where
            qpc_col(icol,:) = 0
            qps_col(icol,:) = 0
            qpg_col(icol,:) = 0
         end where
      end do ! loop over columns
   end subroutine


   subroutine gen_subcol_prec_scops(precip_frac,qpc_avg,qps_avg,qpg_avg,qc_col,qpc_col,qps_col,qpg_col,precip_mode)
      real,intent(in) :: precip_frac(:) ! precipitation fraction
      real,intent(in) :: qpc_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qps_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qpg_avg(:) ! average precipitation mixing ratio
      real,intent(in) :: qc_col(:,:) ! column total cloud mixing ratio
      real,intent(inout) :: qpc_col(:,:) ! column precipitation mixing ratio
      real,intent(inout) :: qps_col(:,:) ! column precipitation mixing ratio
      real,intent(inout) :: qpg_col(:,:) ! column precipitation mixing ratio
      integer,intent(in) :: precip_mode 

      real,allocatable :: ls_p_rate(:,:),cv_p_rate(:,:),cloud_flag(:,:,:),precip_flag(:,:,:)
      integer :: ncol,nlev,icol,ilev
      
      ! get dimension sizes
      ncol = size(qc_col,1)
      nlev = size(qc_col,2)
      
      ! allocate variables for use with precip_scops
      allocate( &
         ls_p_rate(1,nlev),cv_p_rate(1,nlev), &
         cloud_flag(1,ncol,nlev),precip_flag(1,ncol,nlev) &
      )

      ! copy precipitation information to precip_scops arrays
      ls_p_rate(1,:) = qpc_avg + qps_avg + qpg_avg
      cv_p_rate = 0.0
      do icol = 1,ncol
         do ilev = 1,nlev
            if (qc_col(icol,ilev)>0) then
               cloud_flag(1,icol,ilev) = 1
            else
               cloud_flag(1,icol,ilev) = 0
            end if
         end do
      end do

      ! call precip_scops to assign precipitating/non-precipitating columns
      call prec_scops( &
         1,nlev,ncol,ls_p_rate(1,nlev:1:-1),cv_p_rate(1,nlev:1:-1), &
         cloud_flag(1,:,nlev:1:-1),precip_flag(1,:,nlev:1:-1) &
      ) ! precip_scops expects inputs from TOA to surface

      ! adjust precipitation
      if (precip_mode == 1) then ! constrain precip by input precip fraction
         call constrain_subcol_precip(precip_frac,precip_flag(1,:,:))
      else if (precip_mode == 2) then ! remove non-cloudy precip
         where(qc_col == 0)
            precip_flag(1,:,:) = 0
         endwhere
      end if

      ! assign precip mixing ratios
      do icol = 1,ncol
         where (precip_flag(1,icol,:) == 1)
            qpc_col(icol,:) = qpc_avg(:)
            qps_col(icol,:) = qps_avg(:)
            qpg_col(icol,:) = qpg_avg(:)
         else where
            qpc_col(icol,:) = 0
            qps_col(icol,:) = 0
            qpg_col(icol,:) = 0
         end where
      end do ! loop over columns

      deallocate(ls_p_rate,cv_p_rate,cloud_flag,precip_flag)
   end subroutine gen_subcol_prec_scops


   subroutine constrain_subcol_precip(precip_frac,precip_flag)
      real,intent(in) :: precip_frac(:) ! precipitation fraction (TOA to surface)
      real,intent(inout) :: precip_flag(:,:) ! precipitation flag (precipitating or not; TOA to surface)
      integer :: ncol,nlev,icol,ilev

      ! get dimension sizes
      ncol = size(precip_flag,1)
      nlev = size(precip_flag,2)
 
      ! walk down the levels (TOA to surface)
      do ilev = 1,nlev
         if (1.0*sum(precip_flag(:,ilev))/ncol > precip_frac(ilev)) then
            ! need to remove precipitating columns
            icol = 1
            do while (1.0*sum(precip_flag(:,ilev))/ncol > precip_frac(ilev) .and. icol<=ncol)
               if (precip_flag(icol,ilev) == 1) then
                  precip_flag(icol,ilev) = 0
               end if
               icol = icol + 1
            end do
         else if (1.0*sum(precip_flag(:,ilev))/ncol < precip_frac(ilev)) then
            ! need to add precipitating columns
            icol = 1
            do while (1.0*sum(precip_flag(:,ilev))/ncol < precip_frac(ilev) .and. icol<=ncol)
               if (precip_flag(icol,ilev) == 0) then
                  precip_flag(icol,ilev) = 1
               end if
               icol = icol + 1
            end do 
         end if
      end do
   end subroutine constrain_subcol_precip

end module subcol_generator
