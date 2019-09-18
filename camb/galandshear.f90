!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%										%%
!%%	Shear auto and cross correlations output as l(l+1)/2pi C_l^{EE},	%%
!%%	i.e., E-mode shear spectra.						%%
!%%										%%
!%%	**************ALL SPECTRA ARE DIMENSIONLESS!!**************		%%
!%%										%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module CSGalCalc

  use Precision

  use ModelParams
  use Transfer
  use SpherBessels
  use lvalues
  use ModelData
  use constants , ONLY : G !CLUSTERS

  implicit none


  private


  ! All distances and time should be in unit of Mpc, not h^-1 Mpc!!!!!!!

  logical :: DoCS, DoGal, DoCSXGal, Using_CosmoMC
  logical,save :: init = .false.

  character(LEN=40) ShearFileName,GalPSFileName,ShearXGalPSFileName
  character(LEN=40) FracOfGalFileName,RedshiftBinFileName
  character(LEN=40) LMinValFileName,LMaxValFileName
  character(LEN=100) setname,datarep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Tomography stuff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Last index: 1 = shear, 2 = galaxy power spectrum

  integer :: nTomoBin(1:2) ! number of shear/galaxy power spectrum tomography bins 

  real(dl) :: zph_low(nTomoBin_max,1:2),zph_high(nTomoBin_max,1:2)   ! lower and upper limit of redshift bins

  real(dl) :: photo_error		! photo-z scatter: sigma_z = photo_error*(1+z)

  real(dl) :: FracOfGal(nTomoBin_max,1:2) ! fraction of total galaxies in each tomography bin

  logical :: auto_binning(1:2) = .false. ! whether to set redshift bins automagically

  real(dl) :: auto_binning_zmin(1:2),auto_binning_zmax(1:2) ! maximum and minimum redshift, if using autobinning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Redshift distribution parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dl) :: nz_alpha_lss,nz_beta_lss,nz_z0_lss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Maximum multipoles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: lmax_lss(1:2)
  integer :: lmaxindex(1:2)
  integer :: lmin_lss(1:2)
  integer :: lminindex(1:2)
  integer,dimension(:),allocatable :: gal_bin_lmax
  integer,dimension(:),allocatable :: shear_bin_lmax 
  integer,dimension(:),allocatable :: gal_bin_lmin
  integer,dimension(:),allocatable :: shear_bin_lmin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Time integration, etc..!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: nTimeSteps=3000		! Number of time steps for time integration
  ! Tweak this to improve accuracy

  real(dl) :: zsteps(nTimeSteps),chisteps(nTimeSteps),dchisteps(nTimeSteps)
  ! Time steps in terms of redshift and chi 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lensing weights !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  integer :: CS_num_z_tf	! Number of transfer functions used in cosmic shear calculations.
  ! To be set in 	subroutine Transfer_SetForNonlinearLensing_CS(P).

  real(dl) :: weights(nTomoBin_max,nTimeSteps,1:2)	! Lensing and power spectrum weights
  !	real(dl) :: cmblens_weights(nTimeSteps)	! CMB lensing weights

  real(dl) :: galNorm(nTomoBin_max,1:2)	! Normalisation for galaxy distribution function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Linear transfer function !!!!!!!!!!!!!!!!!!!!

  integer :: num_q_trans
  real(dl), dimension(:), allocatable :: q_trans
  real(dl), dimension(:,:), allocatable :: LinTransferData
  real(dl), dimension(:,:), allocatable :: LinTransferDataPsi

  Type (ClTransferData) ScalarTrans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Nonlinear Ratio !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  Type (MatterPowerData) NL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dl), dimension(:,:,:,:), allocatable :: source ! Lensing and galaxy power spectrum source terms
  real(dl), dimension(:,:), allocatable :: cmblens_source ! CMB lensing source term

  real(dl), dimension(:,:,:,:), allocatable :: Source_k 	! k-interpolated source term


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! k integration stuff and "transfer function" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type CSTransfer
     integer :: num_k		! number of k steps for integration
     real(dl), dimension(:), pointer :: dksteps,ksteps	! steps for k integration
     real(dl), dimension(:,:,:,:), pointer :: Delta_p_l_k	! =W(chi)*T(k,chi)*j_l(k*chi) 
  end Type CSTransfer

  Type (CSTransfer) CST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Output data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type OutputCls
     real(dl), dimension(:,:,:,:), allocatable :: CS_Cls
     real(dl), dimension(:,:,:,:), allocatable :: Gal_Cls
     real(dl), dimension(:,:,:,:), allocatable :: CSXGal_Cls
  end Type OutputCls

  Type (OutputCls) :: dummy_cls

  real(dl), dimension(:,:,:,:), allocatable :: iCl_lensing	! Lensing Cl's, at selected ls

  real(dl), dimension(:,:,:,:), allocatable :: Cl_lensing_out	! interpolated Cls

!!!!!!!!!!!!!!!!!!!!!!!!!!! Some CMB stuff for SimData files !!!!!!!!!!!!!!!!!!!!!!!!
!  real(dl), dimension(:), allocatable :: Sigma_T,Sigma_P,fwhm_arcmin

  public Transfer_SetForNonlinearLensingCS,ShePowFoCDriver,SetTomoBins, &
       Cl_lensing_out,GetShePowFoCCls,ShePowFoCInitial, &
       ShearFileName,GalPSFileName,ShearXGalPSFileName,FracOfGalFileName, &
       RedshiftBinFileName,OutputCls,nz_alpha_lss,nz_beta_lss,nz_z0_lss,auto_binning, &
       auto_binning_zmin,auto_binning_zmax,LMinValFileName,LMaxValFileName, &
       dummy_cls,CS_num_z_tf,num_q_trans,q_trans,LinTransferData, &
       GetMatterTransfer,FreeMatterTransferMemory,setname,datarep,OutputSimDataCMB
!       Sigma_T,Sigma_P,fwhm_arcmin,OutputSimDataCMB


contains


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Set tomography bins.						%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

  subroutine SetTomoBins(P)
    implicit none
    type(CAMBParams) :: P     
    integer :: i,j

    nTomoBin = P%nTomoBin
    photo_error = P%photo_error

    nz_z0_lss = P%nz_z0_lss
    nz_alpha_lss = P%nz_alpha_lss
    nz_beta_lss = P%nz_beta_lss

    if (.not. auto_binning(1)) then
       zph_low(1:maxval(P%nTomoBin(:)),1) = P%zph_low(1:maxval(P%nTomoBin(:)),1)
       zph_high(1:maxval(P%nTomoBin(:)),1) = P%zph_high(1:maxval(P%nTomoBin(:)),1)
    else
       call get_redshift_bins(nTomoBin(1),1)
    end if

    if (.not. auto_binning(2)) then
       zph_low(1:maxval(P%nTomoBin(:)),2) = P%zph_low(1:maxval(P%nTomoBin(:)),2)
       zph_high(1:maxval(P%nTomoBin(:)),2) = P%zph_high(1:maxval(P%nTomoBin(:)),2)
    else
       call get_redshift_bins(nTomoBin(2),2)
    end if

  end subroutine SetTomoBins


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Change momentum and redshift sampling of the transfer		%%
  !%%	functions here for better accuracy.				%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! This subroutine sets the redshift and k sampling for the cosmic shear 
  ! calculations.  It replaces the subroutine Transfer_SetForNonlinearLensing(P)
  subroutine Transfer_SetForNonlinearLensingCS(P)
    implicit none
    Type(TransferParams) :: P
    integer :: i, spacing,num_extra_z,istart
    real(dl) :: zmax


    P%kmax = 12.0d0
    P%k_per_logint = 0

    zmax = max(zph_high(nTomoBin(1),1),zph_high(nTomoBin(2),2))
    !if running with clusters only:
    if (zmax .eq. 0.d0)  zmax = 3.d0


    call get_TimeSteps_z(zmax)

    CS_num_z_tf=60
    ! number of redshift values used for cosmic shear calculations

    num_extra_z=nint((9.0-zmax)*AccuracyBoost)
    ! number of extra redshift values, i.e., z>zph_high(ntomobin).
    ! These are used in CAMB for nonlinear corrections in the CMB lensing 
    ! calculations.

    P%num_redshifts=CS_num_z_tf+num_extra_z	! Total number of redshifts


    if (P%num_redshifts > max_transfer_redshifts) &
         stop 'Transfer_SetForNonlinearLensing: Too many redshifts'


    ! Assign extra redshift value
    do i=1,num_extra_z
       P%redshifts(i)=real(num_extra_z-i)/ &
            (num_extra_z/(9.0-zmax))+zmax+1.0
    end do


    ! Assign redshift values for cosmic shear calcs

    spacing=(nTimeSteps-2)/(CS_num_z_tf-2)

    istart=num_extra_z+1

!Luci
    do i=istart,P%num_redshifts-1
       P%redshifts(i) = zsteps(nTimeSteps-(i-istart)*spacing)
    end do

    P%redshifts(P%num_redshifts)=zsteps(1)

  end subroutine Transfer_SetForNonlinearLensingCS


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%								%%%
  !%%	This is the driver!!					%%%
  !%%								%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ShePowFoCDriver(CS,Gal,CSXGal,P,OutCls)
    implicit none

    type (CAMBParams) :: P
    type (OutputCls) :: OutCls
    integer :: i,ik,j,test1,test2
    logical :: CS,Gal,CSXGal
    integer :: time_array_s(8),time_array_e(8)
    real(dl) :: s_time,e_time
    ! All distances and time should be in unit of Mpc, not h^-1 Mpc!!!!!!!
    s_time = time_array_s(5)*3600+time_array_s(6)*60+time_array_s(7)+0.001*time_array_s(8)
    if(FeedbackLevel>0) write(*,*) 'Starting cosmic shear calcs...'

    Using_CosmoMC = P%CSGalCosmoMC

    if(.not. P%flat) then
       write (*,*)  "Sorry, can't handle non-flat models just yet...",P%flat,P%omegak
       pause
    end if

    DoCS = CS
    DoGal = Gal
    DoCSXGal = CSXGal

    if (DoCSXGal .and. .not. (DoCS .and. DoGal)) then
       Print*,'CS/Gal cross-correlation requires calculation of CS and Gal auto-power spectra'
       stop
    end if

    lmax_lss(:) = 0
    lmaxindex(:) = 1
    lmin_lss(:) = 0
    lminindex(:) = 1

    if (DoCS) then
       lmax_lss(1) = P%lmax_CS
       lmin_lss(1) = P%lmin_CS
       do while (lsamp%l(lmaxindex(1)) .lt. lmax_lss(1)) 
          lmaxindex(1) = lmaxindex(1) + 1
       end do
       do while (lsamp%l(lminindex(1)) .lt. lmin_lss(1)) 
          lminindex(1) = lminindex(1) + 1
       end do
    end if


    if (DoGal) then
       lmax_lss(2) = P%lmax_gal
       lmin_lss(2) = P%lmin_gal
       do while (lsamp%l(lmaxindex(2)) .lt. lmax_lss(2)) 
          lmaxindex(2) = lmaxindex(2) + 1
       end do
       do while (lsamp%l(lminindex(2)) .lt. lmin_lss(2)) 
          lminindex(2) = lminindex(2) + 1
       end do
    end if

    if (maxval(lmax_lss(:)) .gt. lsamp%l(lsamp%l0)) stop  'lmax_CS/lmax_gal may not be larger than l_max_scalar'

!!!!!!!!!!!!!!!!!!!!!!! Some initial stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    if(.not. init) then
       call ShePowFoCInitial ! Calculate weights, galaxy distribution, etc..
       init = .true.
    end if

!!!!!!!!!!!!!!!!!!! The real calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call get_TimeSteps_chi  ! Given time steps in terms of z, generate time steps in terms of chi

    call get_weights   ! Calculate lensing/power spectrum weights


    call GetMatterTransfer(P)   ! Get the Matter Transfer functions from CAMB


    call GetShePowFoCCls(P,OutCls)
   
    if(.not. P%CSGalCosmoMC) then

       if(DoCS .and. (ShearFileName .ne. '')) then
          open(unit=34,file=ShearFileName,form='formatted') 
          ! print*, 'COCOS:', ShearFileName
          do i=lmin_lss(1),lmax_lss(1)
             write(34,'(1I5,1000E17.8)') i, OutCls%CS_Cls(i,1:nTomoBin(1),1:nTomoBin(1),1)
           !  print*, i, OutCls%CS_Cls(i,1:nTomoBin(1),1:nTomoBin(1),1)
  
        end do

          close(unit=34)
       end if

       if(DoGal .and. (GalPSFileName .ne. '')) then
          open(unit=35,file=GalPSFileName,form='formatted') 

          do i=lmin_lss(2),lmax_lss(2)
             write(35,'(1I5,1000E17.8)') i, OutCls%Gal_Cls(i,1:nTomoBin(2),1:nTomoBin(2),1)
          end do

          close(unit=35)
       end if

       if(DoCSXGal .and. (ShearXGalPSFileName .ne. '')) then
          open(unit=36,file=ShearXGalPSFileName,form='formatted') 

          do i=maxval(lmin_lss(:)),minval(lmax_lss(:))
             write(36,'(1I5,1000E17.8)') i, OutCls%CSXGal_Cls(i,1:nTomoBin(1),1:nTomoBin(2),1)
          end do

          close(unit=36)
       end if

       if(FracOfGalFileName .ne. '') then
          open(unit=23,file=FracOfGalFileName, form='formatted')          
          write(23,'(1000E17.8)') FracOfGal(1:nTomoBin(1),1)
          write(23,'(1000E17.8)') FracOfGal(1:nTomoBin(2),2)
          close(unit=23)
       end if

       if(RedshiftBinFileName .ne. '') then
          open(unit=67,file=RedshiftBinFileName, form='formatted')


          if (DoCS) then
             write(67,'(100E17.8)') zph_low(1:nTomoBin(1),1)

             write(67,'(100E17.8)') zph_high(1:nTomoBin(1),1)
          else
             write(67,*) ' '
             write(67,*) ' '
          end if
          if (DoGal) then
             write(67,'(100E17.8)') zph_low(1:nTomoBin(2),2)
             write(67,'(100E17.8)') zph_high(1:nTomoBin(2),2)
          end if

          close(unit=67)
       end if

       if(LMinValFileName .ne. '') then
          open(unit=68,file=LMinValFileName, form='formatted')

          allocate(shear_bin_lmin(1:nTomoBin(1)))
          shear_bin_lmin = lmin_lss(1)

          if (DoCS) then
             write(68,'(100I5)') shear_bin_lmin(1:nTomoBin(1))
          else
             write(68,*) ' ' 
          end if

          allocate(gal_bin_lmin(1:nTomoBin(2)))
          gal_bin_lmin = lmin_lss(2)

          if (DoGal) then
             write(68,'(100I5)') gal_bin_lmin(1:nTomoBin(2))
          end if

          close(unit=68)
       end if


       if(LMaxValFileName .ne. '') then
          open(unit=69,file=LMaxValFileName, form='formatted')

          if (DoCS) then
             allocate(shear_bin_lmax(1:nTomoBin(1)))
             shear_bin_lmax = lmax_lss(1)
             write(69,'(100I5)') shear_bin_lmax(1:nTomoBin(1))
          else
             write(69,*) ' ' 
          end if

          if (DoGal) then
             write(69,'(100I5)') gal_bin_lmax(1:nTomoBin(2))
          end if

          close(unit=69)
       end if

!!!!!!!!!!!! Output SimData files if requested !!!!!!!!!!!!
        if(.not. P%CSGalCosmoMC .and. P%OutputSimDataFiles) then

                if (DoCS) then
                        call OutputSimDataCS(P,OutCls)
                end if

                if (DoGal) then
                        call OutputSimDataGal(P,OutCls)
                end if  

                if (DoCSXGal) then
                        call OutputSimDataCSXGal(P,OutCls)
                end if

        end if


    end if


!!!!!!!!!!!!!! Free memory !!!!!!!!!!!!!!!!!!
    call FreeMatterTransferMemory


  end subroutine ShePowFoCDriver








  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
  !%%										%%	
  !%%	Initial stuff	      				 			%%
  !%%										%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ShePowFoCInitial

    implicit none

    if (DoCS) call normalise_galaxy(1)	 ! normalise galaxy distribution functions
    if (DoGal) call normalise_galaxy(2) ! normalise galaxy distribution functions

  end subroutine ShePowFoCInitial

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%										%%
  !%%	The following functions define the galaxy redshift distribution.	%%
  !%% 	Alter them according the survey at hand.				%%
  !%%										%%
  !%%	1. galaxy_final(i,z)							%%
  !%%	2. galaxy_binned_unnorm(i,z)						%%
  !%%	3. galaxy_basic(z)							%%
  !%%	4. normalise_galaxy                                                     %%
  !%%   5. getbins                                                              %%
  !%%										%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function galaxy_final(i,z,n)
    implicit none
    integer :: i,n
    real(dl) :: z,galaxy_final

    galaxy_final = galaxy_binned_unnorm(i,z,n)/galNorm(i,n)

  end function galaxy_final


  !========================================================================

  ! Unnormalised galaxy distribution, AFTER binning
  function galaxy_binned_unnorm(i,z,n)
    implicit none
    integer :: i,n 	
    real(dl) :: z,galaxy_binned_unnorm
    real(dl), parameter :: SQRTTWO=1.4142135623730951d0
    real(dl) :: xhigh,xlow,sig,ngal,blah

    if(photo_error .eq. 0.0d0) then

       ! if sharp edges 
       if((z .ge. zph_low(i,n)) .and. (z .le. zph_high(i,n))) then
          blah=galaxy_basic(z)
       else 
          blah=0.0d0
       endif

    else

       ! photo-z binning a la Ma, Hu & Huterer
       sig = SQRTTWO*photo_error*(1.d0+z)

       xhigh = (zph_high(i,n)-z)/sig
       xlow = (zph_low(i,n)-z)/sig

       blah = 0.5d0*galaxy_basic(z)*(erf(xhigh)-erf(xlow))

    end if

    galaxy_binned_unnorm = blah

  end function galaxy_binned_unnorm

  !=======================================================================

  ! Basic galaxy distribution, BEFORE binning and normalisation.
  function galaxy_basic(z)
    implicit none
    real(dl) :: z,galaxy_basic
    real(dl) :: x

    x=z/nz_z0_lss

    galaxy_basic = (x**nz_alpha_lss)*exp(-x**nz_beta_lss)

  end function galaxy_basic


  !========================================================================
  !z*n(z), needed for calculating average redshift of galaxies per bin
  function z_galaxy_basic(z)
    implicit none
    real(dl) :: z,z_galaxy_basic
    real(dl) :: x

    x=z/nz_z0_lss

    z_galaxy_basic = z*(x**nz_alpha_lss)*exp(-x**nz_beta_lss)

  end function z_galaxy_basic

  !========================================================================

  ! normalise galaxy distribution and calculate the fraction of total
  ! galaxies in each bin 
  subroutine normalise_galaxy(n)
    implicit none
    integer :: i,n	
    real(dl) :: sum

    sum = 0.0d0

    do i=1,nTomoBin(n)
       if(photo_error .eq. 0.d0) then
          ! sharp edges
          galNorm(i,n)=rombint_mod(galaxy_binned_unnorm,i,n,zph_low(i,n),zph_high(i,n),1.0d-6)
       else 
          galNorm(i,n)=rombint_mod(galaxy_binned_unnorm,i,n,zph_low(1,n),zph_high(nTomoBin(n),n),1.0d-6)
       end if

       sum = sum + galNorm(i,n)
    end do

    do i=1,nTomoBin(n)
       FracOfGal(i,n) = galNorm(i,n)/sum
    end do

  end subroutine normalise_galaxy


  !========================================================================


  ! specify galaxy bias as a function of redshift here
  function bias(z)
    implicit none
    real(dl) :: bias
    real(dl), intent(in) :: z

    bias = 1.d0
    !          bias = 1.d0*(1.d0 + z)

  end function bias


  !========================================================================


  subroutine get_redshift_bins(numbins,n)
    implicit none
    integer :: n,numbins
    integer :: i
    real(dl) :: p,ztmp


    zph_low(1,n) = auto_binning_zmin(n)

    do i=1,numbins-1
       p = dble(i)/numbins
       ztmp = quantile(p,galaxy_basic,auto_binning_zmin(n),auto_binning_zmax(n))
       zph_high(i,n) = ztmp
       zph_low(i+1,n) = ztmp
    end do

    zph_high(numbins,n) = auto_binning_zmax(n)

  end subroutine get_redshift_bins





  !========================================================================

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Time Steps for time integration.				%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  ! Time Steps in terms of z
  subroutine get_TimeSteps_z(zmax)
    implicit none
    integer :: i
    real(dl) :: binwidth,ahigh,alow,a,zmax

    alow = 1.0d0
    ahigh = 1.0d0/(1.0d0+zmax)

    binwidth = (log(ahigh)-log(alow))/real(nTimeSteps-1)

    zsteps(1) = 0.0d0
    chisteps(1) = 0.0d0

    do i=2,nTimeSteps
       a = alow*exp(real(i-1)*binwidth)

       zsteps(i) = 1.0d0/a - 1.0d0
    end do

  end subroutine get_TimeSteps_z


  !==========================================================================

  ! Time Steps in terms of chi
  ! All distances and time should be in unit of Mpc, not h^-1 Mpc!
  subroutine get_TimeSteps_chi
    implicit none
    integer :: i

    chisteps(1) = 0.0d0

    do i=2,nTimeSteps
       chisteps(i)=-DeltaTime(1._dl,1._dl/(zsteps(i)+1._dl)) 
    end do

    dchisteps(1)=0.5d0*(chisteps(2)-chisteps(1))

    do i=2,nTimeSteps-1	
       dchisteps(i)=0.5d0*(chisteps(i+1)-chisteps(i-1))
    end do

    dchisteps(nTimeSteps)=0.5d0*(chisteps(nTimeSteps)-chisteps(nTimeSteps-1))

  end subroutine get_TimeSteps_chi


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Lensing weights							%%
  !%%	Also CMB lensing weights.					%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Get lensing/power spectrum weights
  subroutine get_weights
    implicit none
    integer :: i,j,k
    real(dl) :: yout,a,dtauda,w_integrand

    weights(:,:,:) = 0.0d0

    if (DoCS) then

       !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
       !$OMP & PRIVATE(j,yout,w_integrand,a,i,k)

       do j=2,nTimeSteps

          do k=j,nTimeSteps

             a=1.0d0/(1.0d0+zsteps(k))

             w_integrand = fK(chisteps(k)-chisteps(j))/(fK(chisteps(k))*fK(chisteps(j)))
             w_integrand = w_integrand/(dtauda(a)*a*a)

             do i=1,nTomoBin(1)
                yout = w_integrand*galaxy_final(i,zsteps(k),1)
                yout = yout*dchisteps(k)
                weights(i,j,1) = yout + weights(i,j,1)
             end do
          end do

       end do

       !$OMP END PARAllEl DO

    end if

    if (DoGal) then

       !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(STATIC,4) &
       !$OMP & PRIVATE(i,k,a)

       do i=1,nTomoBin(2)
          do k=2,nTimeSteps

             a=1.0d0/(1.0d0+zsteps(k))

             weights(i,k,2) = bias(zsteps(k))*galaxy_final(i,zsteps(k),2)/(dtauda(a)*a*a)

          end do
       end do

       !$OMP END PARAllEl DO

    end if

  end subroutine get_weights


  !==========================================================================

  ! This subroutine outputs f_K(chi), where chi is in units of Mpc, NOT h^-1 Mpc!!	
  ! rofchi(x) is a subroutine from CAMB (in module ModelParams)
  ! CP%r is the factor 1/sqrt(|CP%curv|), where curv = -omegak/(c/h0)^2
  function fK(chi)
    implicit none
    real(dl) :: chi, fK

    fK=rofChi(chi/CP%r)*CP%r

  end function fK

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
  !%%										%%
  !%%	Get the matter transfer functions from CAMB.				%%
  !%%	Transfer functions are stored in MT in CAMB.  Here, we pass the data	%%
  !%%	in MT to variables local to this CS module, so that the data is		%%
  !%%	here in addition.							%%
  !%%										%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine GetMatterTransfer(P)
    implicit none
    integer ijk
    type (CAMBParams) :: P

    num_q_trans = MT%num_q_trans

    print*,num_q_trans,P%Transfer%num_redshifts

    if(allocated(q_trans)) deallocate(q_trans)
    allocate(q_trans(num_q_trans))

    if(allocated(LinTransferData)) deallocate(LinTransferData)
    allocate(LinTransferData(num_q_trans,P%Transfer%num_redshifts))

    if(allocated(LinTransferDataPsi)) deallocate(LinTransferDataPsi)
    allocate(LinTransferDataPsi(num_q_trans,P%Transfer%num_redshifts))

    q_trans = MT%TransferData(1,1:num_q_trans,1)

    LinTransferData(1:num_q_trans,1:P%Transfer%num_redshifts) = &
         MT%TransferData(7,1:num_q_trans,1:P%Transfer%num_redshifts)

!aici 8-->14
    LinTransferDataPsi(1:num_q_trans,1:P%Transfer%num_redshifts) = &
         -MT%TransferData(10,1:num_q_trans,1:P%Transfer%num_redshifts)

    do ijk=1,66
    !print*,P%Transfer%redshifts(ijk),MT%TransferData(10,100,ijk)
    Enddo
  end subroutine GetMatterTransfer

  !================================================================================	

  subroutine FreeMatterTransferMemory

    deallocate(q_trans)
    deallocate(LinTransferData)
    deallocate(LinTransferDataPsi)

  end subroutine FreeMatterTransferMemory

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	The following routines are used to compute the Cls.		%%
  !%%	Uses Limber Approximation at l> 300.				%%
  !%%									%%
  !%%	1. GetShePowFoCCls(P)						%%
  !%%	2. get_source							%%
  !%%	3. get_int_ksteps						%%
  !%%	4. kInterpolateSource						%%
  !%%	5. DoCSFlatIntegration(P,ik)					%%
  !%%	6. locate(x)							%%
  !%%   7. GetClsfromTransfer		                         	%%
  !%%	8. CalcCls					        	%%
  !%%	9. GetCrossClsFromTransfer					%%	
  !%%	10. CalcCrossCls						%%
  !%%	11. CalcCrossCls_Alt						%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine GetShePowFoCCls(P,OutCls)
    implicit none
    type (CAMBParams) :: P
    type (OutputCls) :: OutCls
    integer :: i,nmax,k

    nmax = maxval(nTomoBin)

    if (DoCS .and. (.not. allocated(OutCls%CS_Cls))) then
       !        allocate(OutCls%CS_Cls(2:lSamp%L(lSamp%l0),nTomoBin(1),nTomoBin(1),P%InitPower%nn))
       allocate(OutCls%CS_Cls(2:lmax_lss(1),nTomoBin(1),nTomoBin(1),P%InitPower%nn))
    end if
    if (DoGal .and. (.not. allocated(OutCls%Gal_Cls))) then
       !        allocate(OutCls%Gal_Cls(2:lSamp%L(lSamp%l0),nTomoBin(2),nTomoBin(2),P%InitPower%nn))
       allocate(OutCls%Gal_Cls(2:lmax_lss(2),nTomoBin(2),nTomoBin(2),P%InitPower%nn))
    end if
    if (DoCSXGal .and. (.not. allocated(OutCls%CSXGal_Cls))) then
       !        allocate(OutCls%CSXGal_Cls(2:lSamp%L(lSamp%l0),nTomoBin(1),nTomoBin(2),P%InitPower%nn))
       allocate(OutCls%CSXGal_Cls(2:min(lmax_lss(1),lmax_lss(2)),nTomoBin(1),nTomoBin(2),P%InitPower%nn))
    end if


!!!!!!!!!!!! calculate source terms !!!!!!!!!!!!!!!!!!!!!!

    if (.not. allocated(source)) allocate(source(nmax,num_q_trans,nTimeSteps,1:2))

    call get_source(P)

!!!!!!!!!!!!!!! Get k steps for k integration!!!!!!!!!!!!!!!!!!!
    call get_int_ksteps(P)

!!!!!!!!!!!!!!!! k-interpolate source terms!!!!!!!!!!!!!!!!!!!!!!
    if(allocated(Source_k)) deallocate(Source_k)
    allocate(Source_k(nmax,CST%num_k,nTimeSteps,1:2))

    if(DoCS .or. DoCSXGal) call kInterpolateSource(P,1)
    if(DoGal .or. DoCSXGal) call kInterpolateSource(P,2)

    deallocate(source)

!!!!!!!!!!!!!!!!! Integrate source and bessel function!!!!!!!!!!!!!!!!!!!

    if (associated(CST%Delta_p_l_k)) deallocate(CST%Delta_p_l_k)
    allocate(CST%Delta_p_l_k(nmax,lSamp%l0,CST%num_k,1:2))

    CST%Delta_p_l_k(1:nmax,1:lSamp%l0,1:CST%num_k,1:2) = 0.0d0	! initialisation

    if (DoCS .or. DoCSXGal) then

       !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
       !$OMP & PRIVATE(i)

       do i=1,CST%num_k
          call DoCSFlatIntegration(P,i,1)
       end do

       !$OMP END PARAllEl DO

    end if

    if (DoGal .or. DoCSXGal) then

       !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
       !$OMP & PRIVATE(i)

       do i=1,CST%num_k
          call DoCSFlatIntegration(P,i,2)
       end do

       !$OMP END PARAllEl DO

    end if

    deallocate(Source_k)

!!!!!!!!!!!!!!!Calculate lensing Cls and interpolate !!!!!!!!!!!!!!!!!!!!

    if (DoCS) then
       call GetClsfromTransfer(CST,P,1,1)
       OutCls%CS_Cls(:,:,:,:) = Cl_lensing_out(:,:,:,:)
      ! print*,'buba:'
       !print*, CL_lensing_out(:,:,:,:)
       !print*, '*********************'
    end if

    if (DoGal) then
       call GetClsfromTransfer(CST,P,2,2)
       OutCls%Gal_Cls(:,:,:,:) = Cl_lensing_out(:,:,:,:)
    end if

    if (DoCSXGal) then
       call GetClsfromTransfer(CST,P,1,2)
       OutCls%CSXGal_Cls(:,:,:,:) = Cl_lensing_out(:,:,:,:)
    end if
    Print*, 'Failing'
    deallocate(CST%ksteps)
    deallocate(CST%dksteps)
    deallocate(CST%Delta_p_l_k)

  end subroutine GetShePowFoCCls

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Source terms							%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Basic source term.  Time-interpolated
  subroutine get_source(P)
    implicit none
    Type (CAMBParams) P
    integer :: ibin,ik,ichi,zhi,zlo,z
    real(dl) :: nonlin,transf,a,b,h,rhoMatter,scalefactor
    real(dl) :: transf_Psi

    call GetNLRatio(P) ! Nonlinear ratio
    ! This term here is simply 4 pi G Omega_m rho_cr, i.e., (3/2) Omega_m H_0^2
    rhoMatter=0.5d0*(P%omegac+P%omegab+P%omegan)*grhom


    !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(dynamic) &
    !$OMP & PRIVATE(ichi,ibin,ik,zhi,zlo,z,h,a,b,transf,scalefactor)

    do ichi=1,nTimeSteps

       scalefactor=1.0d0/(1.0d0+zsteps(ichi))

       zlo=1		! These are used for linear interpolation
       zhi=CS_num_z_tf

       do ik=1,num_q_trans	! linear interpolation (in time)

          do while (zhi-zlo .gt. 1)
             z=(zhi+zlo)/2

             if(P%Transfer%redshifts(z) .lt. zsteps(ichi)) then
                zhi=z
             else
                zlo=z
             endif
          end do

          h=P%Transfer%redshifts(zhi)-P%Transfer%redshifts(zlo)

          a=(P%Transfer%redshifts(zhi)-zsteps(ichi))/h
          b=(zsteps(ichi)-P%Transfer%redshifts(zlo))/h

          transf=a*LinTransferData(ik,zlo)+b*LinTransferData(ik,zhi)

          transf_Psi = a*LinTransferDataPsi(ik,zlo)+b*LinTransferDataPsi(ik,zhi)

          nonlin=a*NL%nonlin_ratio(ik,zlo)+b*NL%nonlin_ratio(ik,zhi)
          !nonlin=1.0d0

          ! HaloFit's nonlinear correction can give utterly nonsensical results for non-standard 
          ! models or unusual parameter values.  The following line is a quick fix to keep the 
          ! programme from crashing.
          if (.not. (nonlin .ge. 0.d0))  nonlin = 1.d0

          if (DoCS) then
             do ibin=1,nTomoBin(1)
                source(ibin,ik,ichi,1) = -2.0d0*nonlin*transf_Psi*weights(ibin,ichi,1)
             end do
          end if
          if (DoGal) then
             do ibin=1,nTomoBin(2)
                source(ibin,ik,ichi,2) = nonlin*transf*weights(ibin,ichi,2)*q_trans(ik)**2*(P%H0/100.d0)**2
             end do
          end if

       end do

    end do
    !$OMP END PARAllEl DO

    if(DoGal .and. (.not. Using_CosmoMC)) then 
       call get_bin_max_ell(P)
    end if

    call FreeNonLinearMemory(NL)

  end subroutine get_source


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Get Nonlinear ratio for the transfer function.			%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine GetNLRatio(P)
    use NonLinear, only : NonLinear_GetNonLinRatios
    use InitialPower,only : ScalarPower
    implicit none
    Type (CAMBParams) P

    real(dl) :: hubble,transfsq,matterpower,kh
    integer :: i,iz

    NL%Num_k=num_q_trans
    NL%num_z=P%Transfer%num_redshifts

    call GetNonLinearMemory(NL)

    NL%Redshifts(1:NL%num_z)=P%Transfer%redshifts(1:P%Transfer%num_redshifts)

    ! logkh vector
    hubble=P%H0/100.0d0

    do i=1,NL%num_k	
       NL%log_kh(i)=log(q_trans(i))  !! Suggested by Julien
    end do

    do i=1,NL%num_k
       do iz=1,NL%num_z
          kh=exp(NL%log_kh(i))

          transfsq=LinTransferData(i,iz)**2
          matterpower=transfsq*(ScalarPower(dble(kh*hubble),1)*q_trans(i)*hubble**4*pi*twopi)
          !print*, i, iz, transfsq, ScalarPower(dble(kh*hubble),1), q_trans(i), hubble**4

          ! Suggested by Julien

          NL%matpower(i,iz)=log(matterpower)
       end do
    end do


    ! Spline the logkh-logpk
    call MatterPowerdata_getsplines(NL)
    ! Calculate the ratio sqrt(pk_nonlin/pk_lin)
    call NonLinear_GetNonLinRatios(NL)

  end subroutine GetNLRatio

  !========================================================================

  subroutine GetNonLinearMemory(NL)
    use Transfer, only : MatterPowerData
    implicit none
    Type (MatterPowerData) NL

    allocate(NL%nonlin_ratio(NL%num_k,NL%num_z))
    allocate(NL%log_kh(NL%num_k))
    allocate(NL%matpower(NL%num_k,NL%num_z))	
    allocate(NL%ddmat(NL%num_k,NL%num_z))
    allocate(NL%Redshifts(NL%num_z))

  end subroutine GetNonLinearMemory

  !=========================================================================		

  subroutine FreeNonLinearMemory(NL)
    use Transfer, only : MatterPowerData
    implicit none
    Type (MatterPowerData) NL

    deallocate(NL%nonlin_ratio)
    deallocate(NL%log_kh)
    deallocate(NL%matpower)
    deallocate(NL%ddmat)
    deallocate(NL%Redshifts)

  end subroutine FreeNonLinearMemory

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	k steps for k integration					%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ksteps for k integration
  subroutine get_int_ksteps(P)
    implicit none
    type (CAMBParams) P
    integer :: q_ix
    real(dl) :: sqrtkminus,sqrtkplus,sqrtk
    integer :: numk1,numk2,numk3
    real(dl) :: hubble,kstart,kspacing
    real(dl) :: prelim_dlnk1, prelim_dlnk2, prelim_dlnk3
    !          real(dl), parameter :: prelim_dlnk1=2.0d-1,prelim_dlnk2=1.0d-1,prelim_dlnk3=1.0d-2
    !          real(dl), parameter :: prelim_dlnk1=1.0d-2,prelim_dlnk2=2.0d-3,prelim_dlnk3=1.0d-3
    ! tweak prelim_dlnk1..3 to improve accuracy of k integration.

    !     if (DoCS .and. (.not. DoCSXGal) .and. (.not. DoGal)) then  ! Can get away with relatively low sampling for shear
    !        prelim_dlnk1 = 2.0d-1
    !        prelim_dlnk2 = 1.0d-1
    !        prelim_dlnk3 = 1.0d-2
    !     else ! Galaxy PS requires denser sampling. These settings -> ~0.2% accuracy
    prelim_dlnk1 = 1.d-2
    prelim_dlnk2 = 2.d-3
    prelim_dlnk3 = 1.d-3
    !     end if

    hubble=P%h0/100.0d0

    numk1=nint((-4.0d0-log(hubble*q_trans(1)))/prelim_dlnk1)
    numk2=nint((-2.0d0-(-4.0d0))/prelim_dlnk2)
    numk3=nint((log(hubble*q_trans(num_q_trans))-(-2.0d0))/prelim_dlnk3)+1	

    CST%num_k = numk1 + numk2 + numk3

    if(associated(CST%ksteps)) deallocate(CST%ksteps)
    allocate(CST%ksteps(CST%num_k))

    if(associated(CST%dksteps)) deallocate(CST%dksteps)
    allocate(CST%dksteps(CST%num_k)) 

    kstart = hubble*q_trans(1)
    kspacing = (-4.0d0-log(kstart))/real(numk1,dl)

    do q_ix=1,numk1
       CST%ksteps(q_ix)=kstart*exp((q_ix-1)*kspacing)
    end do

    kstart = exp(-4.0d0)
    kspacing = (-2.0d0-(-4.0d0))/real(numk2,dl)

    do q_ix=numk1+1,numk1+numk2
       CST%ksteps(q_ix) = kstart*exp((q_ix-numk1-1)*kspacing)
    end do

    kstart = exp(-2.0d0)
    kspacing = (log(hubble*q_trans(num_q_trans))- &
         (-2.0d0))/real(numk3-1,dl)

    do q_ix=numk1+numk2+1,CST%num_k
       CST%ksteps(q_ix) = kstart*exp((q_ix-numk1-numk2-1)*kspacing)
    end do

    do q_ix=1,CST%num_k
       sqrtk = sqrt(CST%ksteps(q_ix))

       if(q_ix .ne. 1) then 
          sqrtkminus = sqrt(CST%ksteps(q_ix-1))
       else 
          sqrtkminus = sqrt(CST%ksteps(q_ix))
       end if

       if (q_ix .ne. CST%num_k) then
          sqrtkplus = sqrt(CST%ksteps(q_ix+1))
       else 
          sqrtkplus = sqrt(CST%ksteps(q_ix))
       end if

       CST%dksteps(q_ix)=(sqrtkplus-sqrtkminus)*sqrtk
    end do

  end subroutine get_int_ksteps

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Interpolate source terms in k according to ksteps.		%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Interpolate source terms at more k values.
  subroutine kInterpolateSource(P,n)
    implicit none
    Type (CAMBParams) P
    integer :: q_ix,ibin,itime,n
    integer :: klo,khi,k     
    real(dl) :: hubble,a,b,h,loghi,loglo,blah,logstep
    integer, dimension(:,:), allocatable :: kindex
    real(dl), dimension(:,:), allocatable :: ab

    hubble=P%h0/100.0d0

    allocate(kindex(1:CST%num_k,1:2))
    allocate(ab(1:CST%num_k,1:2))

    do q_ix=1,CST%num_k

       klo=1
       khi=num_q_trans

       logstep=log(CST%ksteps(q_ix))

       do while(khi-klo>1)
          k=(khi+klo)/2

          if(log(hubble*q_trans(k)) > logstep) then
             khi=k
          else
             klo=k
          endif
       end do

       kindex(q_ix,1) = khi
       kindex(q_ix,2) = klo

       h = log(hubble*q_trans(khi)) - log(hubble*q_trans(klo))

       if (h==0.0d0) pause 'bad xa input in linear_interpolate'

       ab(q_ix,1) = (log(hubble*q_trans(khi))-logstep)/h
       ab(q_ix,2) = (logstep-log(hubble*q_trans(klo)))/h

    end do

    !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(dynamic) &
    !$OMP & PRIVATE(itime,q_ix,klo,khi,logstep,loghi,loglo,a,b,h,blah)

    do itime=1,nTimeSteps

       do ibin=1,nTomoBin(n)

          !no need to interpolate if all entries are zero
          if ((maxval(source(ibin,:,itime,n)) .ne. 0.d0) .or.(minval(source(ibin,:,itime,n)) .ne. 0.d0)) then
             do q_ix=1,CST%num_k

                if(source(ibin,kindex(q_ix,2),itime,n)>0 .and. source(ibin,kindex(q_ix,1),itime,n)>0) then
                   blah = ab(q_ix,1)*log(source(ibin,kindex(q_ix,2),itime,n)) + &
                        ab(q_ix,2)*log(source(ibin,kindex(q_ix,1),itime,n))
                   Source_k(ibin,q_ix,itime,n) = exp(blah)
                else
                   Source_k(ibin,q_ix,itime,n) = ab(q_ix,1)*source(ibin,kindex(q_ix,2),itime,n) + &
                        ab(q_ix,2)*source(ibin,kindex(q_ix,1),itime,n)
                end if

             end do
          end if

       end do

    end do

    !$OMP END PARAllEl DO

    deallocate(kindex)
    deallocate(ab)

  end subroutine kInterpolateSource

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%										%%
  !%%	Time integration of Source*Bessel function for one k value		%%
  !%%										%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! CP%flat source integration
  ! Adapted from DoFlatIntegration in CAMB
  subroutine DoCSFlatIntegration(P,ik,n)
    implicit none
    Type (CAMBParams) P
    integer :: ik,n     
    integer j
    real(dl) xlim,xlmax1,k,chif,kappa 
    real(dl) tmin, tmax
    real(dl) a2, J_l, aa(nTimeSteps), fac(nTimeSteps)
    real(dl) xf, sums(nTomoBin(n)), cmblens_sum, dchi
    integer bes_ix, m, bes_index(nTimeSteps),ibin,mplus1,i

    ! Find the position in the xx table for the x correponding to each
    ! timestep

    k=CST%ksteps(ik)

    if (k .le. 1.d0) then
       do j=1,nTimeSteps !Precompute arrays for this k
          xf = abs(k*chisteps(j))
          bes_index(j) = Ranges_indexOf(BessRanges,xf)
          !Precomputed values for the interpolation
          bes_ix = bes_index(j)
          fac(j) = BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
          aa(j) = (BessRanges%points(bes_ix+1)-xf)/fac(j)
          fac(j) = fac(j)**2*aa(j)/6
       end do
    end if

    do j=1,lSamp%l0
       ! Source terms.  Do full calculations up unless l>1000 or 
       ! if k>1.0 Mpc^-1.  Otherwise, use Limber approximation. 
       ! CMB lensing source terms calculated here for low l and for 
       ! small k may not be accurate.  But that is OK, because we will only 
       ! be using these source terms at large k.
       ! See descriptions for CalcCrossCls and CalCrossCls_Alt. 
       if(.not. ((lSamp%l(j) .gt. min(lmax_lss(n),1000)) .or. (k .gt. 1.d0))) then
          xlim = xlimfrac*lSamp%l(j)
          xlim = max(xlim,xlimmin)
          xlim = lSamp%l(j)-xlim
          xlmax1 = 100*max(lSamp%l(j),10)*AccuracyBoost
          tmin = P%tau0-xlmax1/k
          tmax = P%tau0-xlim/k
          tmax = min(P%tau0,tmax)
          tmin = max(P%tau0-chisteps(nTimeSteps-1),tmin)

          sums(1:nTomoBin(n)) = 0.0d0                

          do m=max(1,locate(P%tau0-tmax)), min(nTimeSteps,locate(P%tau0-tmin))

             a2 = aa(m)
             bes_ix = bes_index(m) 

             J_l = a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                  *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(m)) !cubic spline

             J_l = J_l*dchisteps(m)

             do ibin=1,nTomoBin(n)
                sums(ibin) = sums(ibin) + Source_k(ibin,ik,m,n)*J_l
             end do

          end do
       else 
          ! Limber approximation.
          chif = lSamp%l(j)/k

          if(chif>chisteps(nTimeSteps)) then
             sums(1:nTomoBin(n)) = 0.0d0
          else 
             m=locate(chif)

             if (m .le. nTimeSteps-1) then
                mplus1 = m+1
                kappa = (chif-chisteps(m))/(chisteps(mplus1)-chisteps(m))
             else
                mplus1 = m
                kappa = 0.5d0
             end if

             do ibin=1,nTomoBin(n)   
                sums(ibin) = (Source_k(ibin,ik,m,n)*(1-kappa)+&
                     kappa*Source_k(ibin,ik,m+1,n))*sqrt(pi/2/lSamp%l(j))/k 
             end do

          end if
       end if

       CST%Delta_p_l_k(1:nTomoBin(n),j,ik,n) = sums(1:nTomoBin(n))
       !print*, CST%Delta_p_l_k(1:nTomoBin(n),j,ik,n)
    end do

  end subroutine DoCSFlatIntegration

  !===============================================================

  ! This routine is adapted from numerical recipes.
  function locate(x)
    implicit none 
    real(dl) :: x
    integer :: locate
    integer :: jl,jm,ju

    jl=0
    ju=nTimeSteps+1

    do while(ju-jl .gt. 1) 
       jm = (ju+jl)/2

       if((chisteps(nTimeSteps) .ge. chisteps(1)).eqv.(x .ge. chisteps(jm))) then
          jl = jm
       else
          ju = jm
       endif
    end do

    if(x .eq. chisteps(1)) then
       locate = 1
    else if(x .eq. chisteps(nTimeSteps)) then
       locate = nTimeSteps-1
    else
       locate = jl
    endif
    return
  end function locate

  ! This routine is adapted from numerical recipes.
  function zlocate(x)
    implicit none 
    real(dl) :: x
    integer :: zlocate
    integer :: jl,jm,ju

    jl=0
    ju=nTimeSteps+1

    do while(ju-jl .gt. 1) 
       jm = (ju+jl)/2

       if((zsteps(nTimeSteps) .ge. zsteps(1)).eqv.(x .ge. zsteps(jm))) then
          jl = jm
       else
          ju = jm
       endif
    end do

    if(x .eq. zsteps(1)) then
       zlocate = 1
    else if(x .eq. zsteps(nTimeSteps)) then
       zlocate = nTimeSteps-1
    else
       zlocate = jl
    endif
    return
  end function zlocate

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Calculate Cls for lensing tomography at selected ls.		%%
  !%%	Results output as C^E_l, where E=E-mode shear			%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine GetClsfromTransfer(CSTrans,P,n1,n2)
    implicit none
    Type (CSTransfer) CSTrans
    Type (CAMBParams) P
    integer :: pix,n1,n2,n12,kpix

!!!!!!!!!!!!!!!!!!Calculate Cls!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    if (allocated(iCl_lensing)) deallocate(iCl_lensing)	
    allocate(iCl_lensing(lSamp%l0,nTomoBin(n1),nTomoBin(n2),P%InitPower%nn))

    call CalcCls(CSTrans,P,n1,n2)

!!!!!!!!!!!!!!!!Interpolate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(allocated(Cl_lensing_out)) deallocate(Cl_lensing_out)
    allocate(Cl_lensing_out(2:lSamp%L(lSamp%l0),nTomoBin(n1),nTomoBin(n2),P%InitPower%nn))

    n12=P%InitPower%nn
    do pix=1,P%InitPower%nn
    kpix=pix
       call InterpolateCSCl(iCl_lensing(:,:,:,pix),Cl_lensing_out(:,:,:,pix),n1,n2,n12,kpix)
    end do

    deallocate(iCl_lensing)

  end subroutine GetClsfromTransfer

  !===========================================================================

  ! This routine calculates the lensing and galaxy Cls.
  subroutine CalcCls(CSTrans,P,n1,n2)
    implicit none
    Type (CSTransfer) CSTrans
    Type (CAMBParams) P

    integer pix,j
    real(dl) apowers
    integer q_ix,ibin,jbin,n1,n2,lmaxi
    real(dl)   ks(CSTrans%num_k),dlnks(CSTrans%num_k),dlnk
    real(dl) lfactor,lfact2

    if (n1 .eq. n2) lmaxi = lmaxindex(n1)
    if (n1 .ne. n2) lmaxi = minval(lmaxindex(:))

    do pix=1,P%InitPower%nn

       do q_ix=1,CSTrans%num_k
          ks(q_ix) = CSTrans%ksteps(q_ix)  
          dlnks(q_ix) = CSTrans%dksteps(q_ix)/ks(q_ix)
       end do

       !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(static,4) &
       !$OMP & PRIVATE(j,q_ix,dlnk,apowers,ibin,jbin,lfactor,lfact2)

       do j=1,lmaxi

          !Integrate dk/k Delta_l_q**2 * Power(k)
          do ibin=1,nTomoBin(n1)
             do jbin=1,nTomoBin(n2)
                iCl_lensing(j,ibin,jbin,pix)=0.0d0	! initialise
             end do
          end do

          do q_ix=1,CSTrans%num_k
             dlnk=dlnks(q_ix)
             apowers=ScalarPower(ks(q_ix),pix)

             if (n1 .eq. n2) then !autocorrelations

                do ibin=1,nTomoBin(n1)
                   do jbin=1,ibin
                      iCl_lensing(j,ibin,jbin,pix) = iCl_lensing(j,ibin,jbin,pix) + &
                           fourpi*apowers*CSTrans%Delta_p_l_k(ibin,j,q_ix,n1) &
                           *CSTrans%Delta_p_l_k(jbin,j,q_ix,n1)*dlnk
                   end do
                end do

             else !cross-correlation

                do ibin=1,nTomoBin(n1)
                   do jbin=1,nTomoBin(n2)
                      iCl_lensing(j,ibin,jbin,pix) = iCl_lensing(j,ibin,jbin,pix) + &
                           fourpi*apowers*CSTrans%Delta_p_l_k(ibin,j,q_ix,n1) &
                           *CSTrans%Delta_p_l_k(jbin,j,q_ix,n2)*dlnk
                   end do
                end do

             end if

          end do

          !autocorrelations are symmetric in (ibin,jbin)

          if (n1 .eq. n2) then
             do ibin=1,nTomoBin(n1)
                do jbin=ibin+1,nTomoBin(n1)
                   iCl_lensing(j,ibin,jbin,pix)=iCl_lensing(j,jbin,ibin,pix)
                end do
             end do
          end if

          !JH: check

          if ((n1 .eq. 1) .or. (n2 .eq. 1))  then

             lfactor = 0.25d0*real(lSamp%l(j)*lSamp%l(j)-1,dl)*real(lSamp%l(j)+2,dl)*real(lSamp%l(j),dl)

             if (n1 .eq. n2) then

                !Convert phi power spectrum to E-mode power spectrum: C_l^E=1/4 l(l^2-1)(l+2) C_l^phi

                do ibin=1,nTomoBin(n1)
                   do jbin=1,nTomoBin(n1)
                      iCl_lensing(j,ibin,jbin,pix) = lfactor*iCl_lensing(j,ibin,jbin,pix)
                   end do
                end do

             else if (n1 .ne. n2) then

                !Convert C_ell^{Phi,delta} to C_ell^{E,delta}

                do ibin=1,nTomoBin(n1)
                   do jbin=1,nTomoBin(n2)
                      iCl_lensing(j,ibin,jbin,pix) = sqrt(lfactor)*iCl_lensing(j,ibin,jbin,pix)
                   end do
                end do

             end if


          end if

       end do
       !$OMP END PARAllEl DO

    end do

  end subroutine CalcCls

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									%%
  !%%	Interpolate Cls in l						%%
  !%%									%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine InterpolateCSCl(iCl_lensing_in,oCl_lensing_out,n1,n2,n12,pix)
    implicit none     
    integer :: pix,n1,n2,lmaxi,n12
    integer il,llo,lhi, xi
    real(dl) ddCl(nTomoBin(n1),nTomoBin(n2),lSamp%l0)
    real(dl) xl(lSamp%l0)
    integer ibin,jbin
    real(dl) a0,b0,ho
    real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
    real(dl) :: iCl_lensing_in(lSamp%l0,nTomoBin(n1),nTomoBin(n2),n12)
    real(dl) :: oCl_lensing_out(2:lSamp%L(lSamp%l0),nTomoBin(n1),nTomoBin(n2),n12)
    


    if (n1 .eq. n2) lmaxi = lmaxindex(n1)
    if (n1 .ne. n2) lmaxi = minval(lmaxindex(:))

    xl=real(lSamp%l(1:lmaxi),dl)

    !Doing interpolation twice for autocorrelation. Slightly wasteful, but not much time spent here anyway.        
    do ibin=1,nTomoBin(n1)
       do jbin=1,nTomoBin(n2)
          call spline(xl,iCl_lensing_in(:,ibin,jbin,pix),lmaxi,cllo,clhi,ddCl(ibin,jbin,:))
       end do
    end do

    llo=1

    do il=2,lSamp%l(lmaxi)
       xi=il

       if ((xi > lSamp%l(llo+1)).and.(llo < lSamp%l(lSamp%l0))) then
          llo=llo+1
       end if

       lhi=llo+1
       ho=lSamp%l(lhi)-lSamp%l(llo)
       a0=(lSamp%l(lhi)-xi)/ho
       b0=(xi-lSamp%l(llo))/ho

       do ibin=1,nTomoBin(n1)               	
          do jbin=1,nTomoBin(n2)
             oCl_lensing_out(il,ibin,jbin,pix)= a0*iCl_lensing_in(llo,ibin,jbin,pix)+ & 
                  b0*iCl_lensing_in(lhi,ibin,jbin,pix)+ &
                  ((a0**3-a0)* ddCl(ibin,jbin,llo) +(b0**3-b0)*ddCl(ibin,jbin,lhi))*ho**2/6              
          end do
       end do

    end do

  end subroutine InterpolateCSCl

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%								%%	
  !%%	Subroutine quantile(p,func,xmin,xmax)			%%
  !%%	calculates the p-quantile of the (not necessarily       %%
  !%%   normalised) positive definite function func on the      %%
  !%%   interval [xmin,xmax], using the bisection method        %%
  !%%								%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function quantile(p,func,xmin,xmax)
    implicit none
    real(dl) :: quantile
    real(dl),intent(in) :: p,xmin,xmax
    real(dl) :: func
    external func
    real(dl) :: norm,fp
    real(dl) :: xacc = 1.d-6

    real(dl) :: xmid,fmid,f,qtmp
    real(dl) :: dx
    integer :: j
    integer,parameter :: maxit = 40

    if ((p .gt. 1.d0) .or. (p .lt. 0.d0)) stop  'Quantile: p must be between zero and one.'

    norm = cdf(func,xmax,xmin)

    fp = p*norm 

    ! this bit is adapted from the Numerical Recipes routine rtbis
    fmid = cdf(func,xmax,xmin) - fp
    f = cdf(func,xmin,xmin) - fp

    if (f*fmid .gt. 0) stop  'fatal error in quantile (rtbis)'

    if (f .lt. 0.d0) then
       qtmp = xmin
       dx = xmax - xmin
    else
       qtmp = xmax
       dx = xmin - xmax
    endif

    j = 0

    do while (.not. (abs(dx).lt.xacc .or. fmid .eq. 0.d0))
       j = j+1
       dx = dx*0.5d0
       xmid = qtmp + dx
       fmid = cdf(func,xmid,xmin) - fp
       if (fmid .le. 0.d0) qtmp = xmid
       if (j .gt. maxit) stop  'quantile: too many bisections in rtbis'
    end do

    quantile = qtmp

  end function quantile

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%								%%
  !%%    Cumulative distribution function of a distribution     %%
  !%%    function f with support  x .ge. xmin, evaluated at x.  %%
  !%%    If f is not normalised, cdf(f) won't be normalised     %% 
  !%%    either!	                                      	%%
  !%%								%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function cdf(f,x,xmin)
    implicit none
    real(dl) :: cdf
    real(dl),intent(in) :: x,xmin
    real(dl) :: f,rombint2
    external f,rombint2
    real(dl) :: tol = 1d-6

    cdf = rombint2(f,xmin,x,tol,20,100)

  end function cdf

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%								%%
  !%%    Get a rough estimate of the multipole where the linear %%
  !%%    galaxy angular power spectrum becomes unreliable  	%%
  !%% 								%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine get_bin_max_ell(P)
    implicit none
    type (CAMBParams) P
    real(dl),dimension(:),allocatable :: zmean,chi_zmean
    integer, dimension(:),allocatable :: zmean_index
    integer :: i,j
    integer :: izlo,izhi,iz,ell
    real(dl) :: rombint2
    external rombint2
    real(dl) :: nz_norm,znz,nonlin_ratio,cv_ratio,h,a,b,hubble
    real(dl) :: tol = 1d-6

    hubble=P%H0/100.0d0

    allocate(gal_bin_lmax(1:nTomoBin(2)))
    allocate(zmean(1:nTomoBin(2)))
    allocate(chi_zmean(1:nTomoBin(2)))
    allocate(zmean_index(1:nTomoBin(2)))

    ! find mean redshift and corresponding chi for each bin
    do i=1,nTomoBin(2)
       nz_norm = rombint2(galaxy_basic,zph_low(i,2),zph_high(i,2),tol,20,100)
       znz = rombint2(z_galaxy_basic,zph_low(i,2),zph_high(i,2),tol,20,100)
       zmean(i) = znz/nz_norm
       zmean_index(i) = zlocate(zmean(i))
       chi_zmean(i) = chisteps(zmean_index(i))
    end do

    ! assuming k*chi = ell, find the ell where the non-linear correction 
    ! reaches 10% of cosmic variance

    ! JH: changed to more optimistic case; find ell where non-linear 
    ! correction is same as cosmic variance.

    do i=1,nTomoBin(2)

       izlo = 1
       izhi = CS_num_z_tf

       do while (izhi-izlo .gt. 1)
          iz=(izhi+izlo)/2

          if(P%Transfer%redshifts(iz) .lt. zmean(i)) then
             izhi=iz
          else
             izlo=iz
          endif
       end do

       h=P%Transfer%redshifts(izhi)-P%Transfer%redshifts(izlo)

       a = (P%Transfer%redshifts(izhi)-zmean(i))/h
       b = (zmean(i)-P%Transfer%redshifts(izlo))/h

       nonlin_ratio = 1.d0
       cv_ratio = 1.d0

       j = 0
       do while (abs(nonlin_ratio-1.d0) .lt. cv_ratio)
          j = j+1   
          ell = nint(hubble*chi_zmean(i)*q_trans(j))
          cv_ratio = sqrt(2.d0/(2*dble(ell)+1.d0))
          nonlin_ratio = a*NL%nonlin_ratio(j,izlo) + b*NL%nonlin_ratio(j,izhi)
       end do

       gal_bin_lmax(i) = min(ell,lmax_lss(2))

    end do


    deallocate(zmean)
    deallocate(chi_zmean)
    deallocate(zmean_index)

  end subroutine get_bin_max_ell

    subroutine OutputSimDataCMB(P)
        implicit none

        type(CAMBParams) :: P

        integer :: nc,ll
        real(dl) :: NoiseCl(2,P%lmin_cmb:P%lmax_cmb)
        real(dl) :: sigma2(P%nchan),fwhm_rad(P%nchan)
        real(dl) :: TT,EE,TE,BB,xlc
        character(LEN=100) :: setoutname,cmboutname,datname

        datname = trim(setname)//"_cmb.dat"
        cmboutname = trim(datarep)//trim(datname)
        setoutname = trim(datarep)//trim(setname)//"_cmb.dataset"
        
        xlc = sqrt(8.d0*log(2.d0))
        fwhm_rad = P%fwhm_arcmin/60.d0*pi/180.d0
        sigma2 = (fwhm_rad/xlc)**2

        NoiseCl = 0.d0
        TT = 0.d0
        TE = 0.d0
        EE = 0.d0
        BB = 0.d0

        do ll = P%lmin_cmb,P%lmax_cmb
                do nc = 1,P%nchan
                        NoiseCl(1,ll) = NoiseCl(1,ll) + &
                                1.d0/((fwhm_rad(nc)*P%Noise_Sigma_T(nc))**2) &
                                *exp(real(-ll*(ll+1),KIND(1.d0))*sigma2(nc))
                        NoiseCl(2,ll) = NoiseCl(2,ll) + &
                                1.d0/((fwhm_rad(nc)*P%Noise_Sigma_P(nc))**2) &
                                *exp(real(-ll*(ll+1),KIND(1.d0))*sigma2(nc))
                end do
        end do

        NoiseCl = 1.d0/NoiseCl

        call CreateTxtFile(cmboutname,34)

        if(P%sim_random_cmb)  stop ('Randomisation of CMB data not implemented yet!')

        do ll = P%lmin_cmb,P%lmax_cmb
           if(.not. P%DoLensing) then
               TT = Cl_scalar(ll,1,C_Temp)
               EE = Cl_scalar(ll,1,C_E)
               TE = Cl_scalar(ll,1,C_Cross)
            else !if doing lensing
               TT = Cl_lensed(ll,1,CT_Temp)
               EE = Cl_lensed(ll,1,CT_E)
               TE = Cl_lensed(ll,1,CT_Cross)
               BB = Cl_lensed(ll,1,CT_B)
            end if

            if(P%WantTensors .and. ll .le. P%Max_l_tensor) then
                TT = TT + Cl_tensor(ll,1,CT_Temp)
                EE = EE + Cl_tensor(ll,1,CT_E)
                TE = TE + Cl_tensor(ll,1,CT_Cross)
                BB = BB + Cl_tensor(ll,1,CT_B)
             end if

             TT = TT*2.d0*pi/real(ll*(ll+1),KIND(1.d0))*(P%tcmb*1.d6)**2 + &
                        NoiseCl(1,ll)
             EE = EE*2.d0*pi/real(ll*(ll+1),KIND(1.d0))*(P%tcmb*1.d6)**2 + &
                        NoiseCl(2,ll)
             TE = TE*2.d0*pi/real(ll*(ll+1),KIND(1.d0))*(P%tcmb*1.d6)**2
             BB = BB*2.d0*pi/real(ll*(ll+1),KIND(1.d0))*(P%tcmb*1.d6)**2 + &
                        NoiseCl(2,ll)

             if (P%sim_ncmbcls .eq. 1) then
                write(34,'(1I5,13E17.8)') ll, TT, NoiseCl(1,ll), P%fskycmb
             else if (P%sim_ncmbcls .eq. 3) then
                write(34,'(1I5,13E17.8)') ll, TT, TE, TE, EE, NoiseCl(1,ll), NoiseCl(2,ll), P%fskycmb
             else if (P%sim_ncmbcls .eq. 4) then
                write(34,'(1I5,13E17.8)') ll, TT, TE, 0.d0, TE, EE, 0.d0, 0.d0, &
                        0.d0, BB, NoiseCl(1,ll), NoiseCl(2,ll), NoiseCl(2,ll), P%fskycmb
             end if
        end do

        call CreateTxtFile(setoutname,37)
        write (37,'(A,A)')  "name = ",trim(setname)
        write (37,'(A,L1)') "all_l_exact = ",.true.
        write (37,'(A,L1)') "CMB = ",.true.
        write (37,'(A,I5)') "sim_ncmbcls = ",P%sim_ncmbcls
        write (37,'(A,I5)') "sim_lmin = ",P%lmin_cmb
        write (37,'(A,I5)') "sim_lmax = ",P%lmax_cmb
        write (37,'(A,A)') "datafile = ",trim(datname)
        write (37,'(A)') "#CAMB settings:" 
        write (37,'(A,L1)') "CMB_lensed = ",P%DoLensing
       	write (37,'(A,I6)') "CAMB_l_max_scalar = ",P%Max_l
        write (37,'(A,E13.5)') "CAMB_k_eta_max_scalar = ",P%Max_eta_k 
        write (37,'(A,I6)') "CAMB_l_max_tensor = ",P%Max_l_tensor
        write (37,'(A,E13.5)') "CAMB_k_eta_max_tensor = ",P%max_eta_k_tensor
        close(37)

    end subroutine OutputSimDataCMB

  subroutine OutputSimDataCS(P,OutCls)
        implicit none
  
        type (CAMBParams) :: P
        type (OutputCls) :: OutCls
 
        real(dl) :: NoiseCov(nTomoBin(1))
        real(dl) :: cov(lmin_lss(1):lmax_lss(1),nTomoBin(1),nTomoBin(1)) 
        integer :: l,i
        character(LEN=100) :: csoutname,setoutname,datname

        call get_noise_spec(P%ngal_lss,P%mean_int_ellip_lss,NoiseCov)

        cov = OutCls%CS_Cls(:,1:nTomoBin(1),1:nTomoBin(1),1)

        do l = lmin_lss(1),lmax_lss(1)
                do i = 1,nTomoBin(1)
                        cov(l,i,i) = cov(l,i,i)+NoiseCov(i) 
                end do
        end do

        datname = trim(setname)//"_shear.dat"
        csoutname = trim(datarep)//trim(datname)

        call CreateTxtFile(csoutname,34)

        do l=lmin_lss(1),lmax_lss(1)
                write (34,'(1I5,1111E17.8)') l,cov(l,1:nTomoBin(1),1:nTomoBin(1)), &
                       NoiseCov(1:nTomoBin(1)),P%fskylss
        end do

        close(34)

        setoutname = trim(datarep)//trim(setname)//"_shear.dataset"

        call CreateTxtFile(setoutname,37)
       
        write (37,'(A,A)')  "name = ",trim(setname)
        write (37,'(A,L1)') "all_l_exact = ",.true.
        write (37,'(A,L1)') "CosmicShear = ",.true.
        write (37,'(A,I5)') "sim_ntomobin_cs = ",nTomoBin(1)
        write (37,'(A,I5)') "sim_lmin = ",lmin_lss(1)
        write (37,'(A,I5)') "sim_lmax = ",lmax_lss(1)
        write (37,'(A,100I5)') "lmin_vector_cs =",shear_bin_lmin(1:nTomoBin(1))
        write (37,'(A,100I5)') "lmax_vector_cs =",shear_bin_lmax(1:nTomoBin(1))
        write (37,'(A,100F12.5)') "zlow_vector_cs =",zph_low(1:nTomoBin(1),1)
        write (37,'(A,100F12.5)') "zhigh_vector_cs =",zph_high(1:nTomoBin(1),1)
        write (37,'(A,F12.5)') "photo_error = ",photo_error
        write (37,'(A,F12.5)') "nz_z0_lss = ",nz_z0_lss
        write (37,'(A,F12.5)') "nz_alpha_lss = ",nz_alpha_lss
        write (37,'(A,F12.5)') "nz_beta_lss = ",nz_beta_lss
        write (37,'(A,A)') "datafile = ",trim(datname)
        close(37)

  end subroutine OutputSimDataCS

  subroutine OutputSimDataGal(P,OutCls)
        implicit none
  
        type (CAMBParams) :: P
        type (OutputCls) :: OutCls
 
        real(dl) :: NoiseCov(nTomoBin(2))
        real(dl) :: cov(lmin_lss(2):lmax_lss(2),nTomoBin(2),nTomoBin(2)) 
        integer :: l,i
        character(LEN=100) :: galoutname,setoutname,datname

        call get_gal_noise_spec(P%ngal_lss,NoiseCov)

        cov = OutCls%Gal_Cls(:,1:nTomoBin(2),1:nTomoBin(2),1)

        do l = lmin_lss(2),lmax_lss(2)
                do i = 1,nTomoBin(2)
                        cov(l,i,i) = cov(l,i,i)+NoiseCov(i) 
                end do
        end do

        datname = trim(setname)//"_gal.dat"
        galoutname = trim(datarep)//trim(datname)

        call CreateTxtFile(galoutname,34)

        do l=lmin_lss(2),lmax_lss(2)
                write (34,'(1I5,1111E17.8)') l,cov(l,1:nTomoBin(2),1:nTomoBin(2)), &
                       NoiseCov(1:nTomoBin(2)),P%fskylss
        end do

        close(34)

        setoutname = trim(datarep)//trim(setname)//"_gal.dataset"

        call CreateTxtFile(setoutname,37)
        
        write (37,'(A,A)')  "name = ",trim(setname)
        write (37,'(A,L1)') "all_l_exact = ",.true.
        write (37,'(A,L1)') "GalaxyPower = ",.true.
        write (37,'(A,I5)') "sim_ntomobin_gal = ",nTomoBin(2)
        write (37,'(A,I5)') "sim_lmin = ",lmin_lss(2)
        write (37,'(A,I5)') "sim_lmax = ",lmax_lss(2)
        write (37,'(A,100I5)') "lmin_vector_gal =",gal_bin_lmin(1:nTomoBin(2))
        write (37,'(A,100I5)') "lmax_vector_gal =",gal_bin_lmax(1:nTomoBin(2))
        write (37,'(A,100F12.5)') "zlow_vector_gal =",zph_low(1:nTomoBin(2),2)
        write (37,'(A,100F12.5)') "zhigh_vector_gal =",zph_high(1:nTomoBin(2),2)
        write (37,'(A,F12.5)') "photo_error = ",photo_error
        write (37,'(A,F12.5)') "nz_z0_lss = ",nz_z0_lss
        write (37,'(A,F12.5)') "nz_alpha_lss = ",nz_alpha_lss
        write (37,'(A,F12.5)') "nz_beta_lss = ",nz_beta_lss
        write (37,'(A,A)') "datafile = ",trim(datname)
        close(37)

  end subroutine OutputSimDataGal

  subroutine OutputSimDataCSXGal(P,OutCls)
        implicit none
  
        type (CAMBParams) :: P
        type (OutputCls) :: OutCls
 
        real(dl) :: cov(maxval(lmin_lss):minval(lmax_lss),nTomoBin(1),nTomoBin(2)) 
        integer :: l,i
        character(LEN=100) :: crossoutname,setoutname,datname,csdatname,galdatname

        cov = OutCls%CSXGal_Cls(:,1:nTomoBin(1),1:nTomoBin(2),1)
        
        csdatname = trim(setname)//"_shear.dat"
        galdatname = trim(setname)//"_gal.dat"
        datname = trim(setname)//"_shearXgal.dat"
        crossoutname = trim(datarep)//trim(datname)

        call CreateTxtFile(crossoutname,34)

        do l=maxval(lmin_lss),minval(lmax_lss)
                write (34,'(1I5,1111E17.8)') l,cov(l,1:nTomoBin(1),1:nTomoBin(2)),P%fskylss
        end do

        close(34)

        setoutname = trim(datarep)//trim(setname)//"_shearXgal.dataset"

        call CreateTxtFile(setoutname,37)
        
        write (37,'(A,A)')  "name = ",trim(setname)
        write (37,'(A,L1)') "all_l_exact = ",.true.
        write (37,'(A,L1)') "CosmicShear = ",.true.
        write (37,'(A,L1)') "GalaxyPower = ",.true.
        write (37,'(A,I5)') "sim_ntomobin_cs = ",nTomoBin(1)
        write (37,'(A,I5)') "sim_ntomobin_gal = ",nTomoBin(2)
        write (37,'(A,I5)') "sim_lmin = ",maxval(lmin_lss)
        write (37,'(A,I5)') "sim_lmax = ",minval(lmax_lss)
        write (37,'(A,100I5)') "lmin_vector_cs =",shear_bin_lmin(1:nTomoBin(1))
        write (37,'(A,100I5)') "lmax_vector_cs =",shear_bin_lmax(1:nTomoBin(1))
        write (37,'(A,100I5)') "lmin_vector_gal =",gal_bin_lmin(1:nTomoBin(2))
        write (37,'(A,100I5)') "lmax_vector_gal =",gal_bin_lmax(1:nTomoBin(2))
        write (37,'(A,100F12.5)') "zlow_vector_cs =",zph_low(1:nTomoBin(1),1)
        write (37,'(A,100F12.5)') "zhigh_vector_cs =",zph_high(1:nTomoBin(1),1)
        write (37,'(A,100F12.5)') "zlow_vector_gal =",zph_low(1:nTomoBin(2),2)
        write (37,'(A,100F12.5)') "zhigh_vector_gal =",zph_high(1:nTomoBin(2),2)
        write (37,'(A,F12.5)') "photo_error = ",photo_error
        write (37,'(A,F12.5)') "nz_z0_lss = ",nz_z0_lss
        write (37,'(A,F12.5)') "nz_alpha_lss = ",nz_alpha_lss
        write (37,'(A,F12.5)') "nz_beta_lss = ",nz_beta_lss
        write (37,'(A,A)') "datafile = ",trim(datname)
        write (37,'(A,A)') "datafile_cs = ",trim(csdatname)
        write (37,'(A,A)') "datafile_gal = ",trim(galdatname)
        close(37)
  end subroutine OutputSimDataCSXGal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%                                                             %%
!%%     Lensing noise power spectrum                            %%
!%%                                                             %%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine get_noise_spec(ngal,int_ellip,NoiseCov)
        implicit none

        !integer, intent(IN) :: nTomoBin(1:2)
        real(KIND(1.d0)), intent(IN) :: ngal ! galaxy surface density, in arcmin^-2
        real(KIND(1.d0)), intent(IN) :: int_ellip ! intrinsic ellipticity.
        !real(KIND(1.d0)), intent(IN) :: FracOfGal(nTomoBin(1),1:2)
        real(KIND(1.d0)), intent(OUT) :: NoiseCov(nTomoBin(1))

        double precision, parameter :: AMIN2STER=8.46159d-8     ! arcmin^2 to steradian conversion factor
        double precision :: noise

        integer :: i

        noise=AMIN2STER/ngal

        do i=1,nTomoBin(1)
                NoiseCov(i)=noise*int_ellip*int_ellip/FracOfGal(i,1)
        end do

        end subroutine get_noise_spec


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%                                                             %%
!%%     Galaxy noise power spectrum                             %%
!%%                                                             %%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine get_gal_noise_spec(ngal,NoiseCov)
        implicit none

        real(KIND(1.d0)), intent(IN) :: ngal ! galaxy surface density, in arcmin^-2
        real(KIND(1.d0)), intent(OUT) :: NoiseCov(nTomoBin(2))

        double precision, parameter :: AMIN2STER=8.46159d-8     ! arcmin^2 to steradian conversion factor
        double precision :: noise

        integer :: i

        noise=AMIN2STER/ngal

        do i=1,nTomoBin(2)
                NoiseCov(i)=noise/FracOfGal(i,2)
        end do

        end subroutine get_gal_noise_spec

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%								  %%	
  !%%	Romberg Integration					  %%
  !%%	This is different from the one in the CAMB package.	  %%
  !%%								  %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !  Rombint returns the integral from a to b of using Romberg integration.
  !  The method converges provided that f(ibin,x) is continuous in (a,b).
  !  f must be real(dl) and must be declared external in the calling
  !  routine.  tol indicates the desired relative accuracy in the integral.
  function rombint_mod(f,ibin,n,a,b,tol)
    use Precision
    implicit none
    integer, parameter :: MAXITER=20,MAXJ=5
    dimension g(MAXJ+1)
    integer :: ibin,n
    real(dl) :: f
    external f
    real(dl) :: rombint_mod
    real(dl), intent(in) :: a,b,tol
    integer :: nint, i, k, jmax, j
    real(dl) :: h, gmax, error, g, g0, g1, fourj

    h=0.5d0*(b-a)
    gmax=h*(f(ibin,a,n)+f(ibin,b,n))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0

10  i=i+1

    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
         go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl

    do k=1,nint
       g0=g0+f(ibin,a+(k+k-1)*h,n)
    end do

    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl

    do j=1,jmax
       !  Use Richardson extrapolation.
       fourj=4._dl*fourj
       g1=g0+(g0-g(j))/(fourj-1._dl)
       g(j)=g0
       g0=g1
    end do

    if (abs(g0).gt.tol) then
       error=1._dl-gmax/g0
    else
       error=gmax
    end if

    gmax=g0
    g(jmax+1)=g0

    go to 10

40  rombint_mod=g0

    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
       write(*,*) 'Warning: Rombint failed to converge; '
       write (*,*)'integral, error, tol:', rombint_mod,error, tol
    end if

  end function rombint_mod

end module CSGalCalc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
