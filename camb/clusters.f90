!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%										%%
!%%	Shear auto and cross correlations output as l(l+1)/2pi C_l^{EE},	%%
!%%	i.e., E-mode shear spectra.						%%
!%%										%%
!%%	**************ALL SPECTRA ARE DIMENSIONLESS!!**************		%%
!%%										%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module Clusters

  use Precision

  use ModelParams
  use Transfer
  use SpherBessels
  use lvalues
  use ModelData
  use constants , ONLY : G !CLUSTERS
  use CSGalCalc

  implicit none

  private

  character(LEN=40) ClRedshiftBinFileName,ClMassBinFileName,ClFileName
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cluster stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: zres = 1000 !redshift sampling before creating bins for fiducial
  integer, parameter :: Mres = 10000 !mass sampling before creating bins for fiducial
  real(dl), parameter :: Mhi = 1.d16 !upper mass limit above which we find essentially no clusters
  real(dl), parameter :: Mlo = 1.d13 !lower mass limit which the detection threshold 'always' exceeds
  real(dl), parameter :: eps_ode = 1.d-6 !precision for solving the differential equations
  real(dl), parameter :: acc_deltami = 1d-8 !accuracy for finding the initial delta corresponding to collapse at z
  real(dl), parameter :: tol_sigmaCR = 1.d-6 !relative tolerance for critical surface mass density
  real(dl), parameter :: tol_alpha = 1.d-5 !relative tolerance for calculation of alpha in eq. (5.2) of arXiv:1304:2321
  real(dl), parameter :: tol_cmf = 1.d-2 !relative tolerance for integration across mass bins
  real(dl), parameter :: acc_Mmin = 1.d10 !accuracy for the detection threshold
  real(dl), parameter :: tol_clustersperbin = 1.d-2 !relative tolerance for integration across redshift bins

  Type OutputClusters
     real(dl), dimension(:,:), allocatable :: clustercount
  end Type OutputClusters

  Type (OutputClusters) :: dummy_clusters

  public ClusterFoCDriver,OutputClusters, &
         ClRedshiftBinFileName,ClMassBinFileName,ClFileName,dummy_clusters


contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routines for the cluster survey                       	  %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine ClusterFoCDriver(P,OutClusters)
    implicit none

    type (CAMBParams) :: P
    type (OutputClusters) :: OutClusters
    real(dl) :: clustersperbin(P%cl_zbins,P%cl_Mbins) !CLUSTERS
    integer :: i,j
    integer :: time_array_s(8),time_array_e(8)
    real(dl) :: s_time,e_time
    if (.not. allocated(OutClusters%clustercount)) then
       allocate(OutClusters%clustercount(P%cl_zbins,P%cl_mbins))
    end if

    call GetMatterTransfer(P)           ! Get the Matter Transfer functions from CAMB

    if (.not. P%CSGalCosmoMC) then
       call GetClusterBinsFiducial(P%cl_RedshiftBins,P%cl_MassBins)
    end if

    call FillClusterBins(clustersperbin,P%cl_RedshiftBins,P%cl_MassBins)
    
    OutClusters%clustercount(1:P%cl_zbins,1:P%cl_mbins) = clustersperbin(1:P%cl_zbins,1:P%cl_mbins)

    if (.not. P%CSGalCosmoMC) then
       
       if (ClFileName .ne. '') then
          open(unit=27,file=ClFileName,form='formatted')
          do i = 1,P%cl_zbins
             write(27,'(1000E17.8)')(clustersperbin(i,j),j=1,P%cl_Mbins)
          end do
          close(unit=27)
       end if

       if (ClRedshiftBinFileName .ne. '') then
          open(unit=28,file=ClRedshiftBinFileName,form='formatted')
          write(28,'(40E17.8)')(P%cl_RedshiftBins(j),j=1,P%cl_zbins+1)
          close(unit=28)
       end if

       if (ClMassBinFileName .ne. '') then
          open(unit=29,file=ClMassBinFileName,form='formatted')
          do i = 1,P%cl_zbins
             write(29,'(50E17.8)')(P%cl_MassBins(i,j),j=1,P%cl_Mbins+1)
          end do
          close(unit=29)
       end if

    end if

    if (.not. P%CSGalCosmoMC .and. P%OutputSimDataFiles) then
        call OutputSimDataCl(P,OutClusters)
    end if

    call FreeMatterTransferMemory

  end subroutine ClusterFoCDriver

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routines for setting up redshift and mass                 %%
  !%%           with approx. equal numbers of clusters in                 %%
  !%%           each bin.                                                 %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !Routine for setting redshift and mass bins such that each contains approx same number of clusters in fiducial model
  subroutine GetClusterBinsFiducial(zb,Mb)
    implicit none
    real(dl) :: zb(cl_zbins_max+1),Mb(cl_zbins_max,cl_Mbins_max+1),factor(CP%cl_zbins+1),pop(zres,Mres)
    integer :: kb(CP%cl_zbins+1)
    real(dl) :: Nb,temp,fac,dz,dlogM,w,Ntot,zlo,zhi,NumClBin
    integer :: i,j,k,l

    !Initialization of parameters - kb is used to keep track of which redshift value bursts the bin and factor is used for linear interpolation between the two z-values that brackets the upper bin limit.
    call Getzlo(zb(1),3)
    kb(1) = 1
    factor(1) = 1.d0
    dz = (CP%zmax_cl-zb(1))/(zres-1.d0)
    dlogM = (log(Mhi)-log(Mlo))/(Mres-1.d0)

    !Call of clusters subroutine which returns the number of clusters on a fine redshift and mass grid - zhi is the lowest redshift at which the detection threshold exceeds Mhi (10^16 Msun).
    call populate(pop,dz,dlogM,zb(1),zhi)

    if (zhi .eq. CP%zmin_cl) stop 'Bad choice of fiducial model causing spherical collapse to fail at all redshifts'

    !The total number of clusters is computed
    Ntot = sum(pop)

    !Loop for distributing redshift bins
    k = 0
    Nb = 0.d0
    do j = 1,CP%cl_zbins-1
       do 
          k = k + 1

          !Summation along mass dimension of the grid
          temp = sum(pop(k,:))
          Nb = Nb + temp

          !Each redshift bin should contain approx. Ntot/#z-bins clusters - exit loop when this is exceeded.
          if (Nb .gt. Ntot/CP%cl_zbins) goto 98
       end do

       !Compute excess of clusters in the redshift bin and preform linear interpolation to determine the upper bin limit
98     Nb = Nb-Ntot/CP%cl_zbins
       factor(j+1) = Nb/temp
       kb(j+1) = k
       zb(j+1) = zb(1)+(k-factor(j+1))*dz
    end do

    !Fill in information about upper limit of the last redshift bin
    zb(CP%cl_zbins+1) = zhi
    factor(CP%cl_zbins+1) = 0.d0
    kb(CP%cl_zbins+1) = zres

    !For every redshift bin - mass bins are determined in the same manner
    do j = 1,CP%cl_zbins
       Mb(j,1) = Mlo
       Nb = 0.d0
       i = 0

       !Loop over mass bins
       do l = 1,CP%cl_Mbins-1
          do
             i = i + 1
             temp = 0.d0
             do k = kb(j),kb(j+1)
                w = 1.d0

                !Make use of the linear interpolation used to determine upper redshift bin limit
                if (k .eq. kb(j)) w = factor(j)
                if (k .eq. kb(j+1)) w = 1-factor(j+1)
                temp = temp + w*pop(k,i)
             end do
             Nb = Nb+temp

             !Each mass bin should contain approx. Ntot/#z-bins/#M-bins clusters
             if (Nb .gt. Ntot/CP%cl_zbins/CP%cl_Mbins) goto 99
          end do

          !Compute excess of clusters and do linear interpolation
99        Nb = Nb-Ntot/CP%cl_zbins/CP%cl_Mbins
          fac = Nb/temp
          Mb(j,l+1) = exp(log(Mb(j,1)) + (i-fac)*dlogM)
       end do

       !Upper limit of the last mass bin
       Mb(j,CP%cl_Mbins+1) = Mhi
    end do

  end subroutine GetClusterBinsFiducial

  !Subroutine for finding the lowest redshift where the detection threshold is above Mhi
  subroutine Getzlo(zlo,factor)
    implicit none
    real(dl) :: zlo,z,delta_c,Delta_vir,Mmin
    integer :: factor,i

    do i = 1,factor*zres-1
       z = CP%zmin_cl+i*(CP%zmax_cl-CP%zmin_cl)/(factor*zres-1.d0)
       
       call spherical_collapse(z,delta_c,Delta_vir)
       
       if (delta_c .gt. 0.d0 .and. Delta_vir .gt. 0.d0) then
          call rtsafe(CC_Mmin,z,Mlo,Mhi,1.d8,Mmin,Delta_vir)
       end if

       if (Mmin .gt. 0.d0) then
          zlo = z
          return
       end if
    end do

    stop 'Error from cluster subroutine: detection threshold outside range of Mlo to Mhi for all redshifts'
  end subroutine Getzlo

  !Function for filling fine equidistant redshift and mass grid with clusters - this is for later division into bins of equal numbers
  subroutine populate(pop,dz,dlogM,zlo,zhi)
    implicit none
    real(dl) :: pop(zres,Mres),Mint(Mres),sigma(Mres),dlnsdlnM(Mres)
    real(dl) :: dz,dlogM,z,Mmin,dV,w,cmf,hmf,zlo,zhi,delta_c,Delta_vir,nbg,rombint
    integer :: i,j,k
    external :: rombint
    !Initialize mass grid
    do j = 1,Mres
       Mint(j) = exp(log(Mlo) + (j-1)*dlogM)
    end do

    !Highest possible redshift - is changed along the way to the lowest redshift at which the detection threshold exceeds Mhi (10^16 Msun)
    zhi = CP%zmax_cl

    !$OMP PARAllEl DO DEFAUlT(SHARED) &
    !$OMP & PRIVATE(i,j,k,delta_c,Delta_vir,z,Mmin,sigma,dlnsdlnM,cmf,w,hmf,dV)

    !Loop over redshift
    do i = 1,zres
       z = zlo+(i-1.d0)*dz
       delta_c = 0.d0
       Delta_vir = 0.d0

       !Solve the spherical collapse to obtain delta_c(z) and Delta_vir(z)
       if (z .ne. 0.) call spherical_collapse(z,delta_c,Delta_vir)

       if (delta_c .gt. 0.d0 .and. Delta_vir .gt. 0.d0 .and. z .ne. 0.) then

          !Based on Delta_vir(z) compute detection threshold
          call rtsafe(CC_Mmin,z,Mlo,Mhi,1.d8,Mmin,Delta_vir) 

          if (Mmin .lt. 0) then
             call rtsafe(CC_Mmin,z,1.d0,Mhi,1.d8,Mmin,Delta_vir)
             if (Mmin .eq. -1.0) stop 'Error from cluster subroutine: Mmin is lower than Mlo - change Mlo in clusters.f90'
          end if

          !Comoving volume element
          dV = comvol(z)

          !Only proceed if detection threshold is smaller than Mhi (10^16 Msun)
          if (Mmin.gt.0.d0) then

             !Compute sigma(M,z) for entire mass grid
             do j = 1,Mres
                sigma(j) = Get_sigma(Mint(j),z)
             end do

             !Calculate the derivative of ln(sigma) with respect to ln(M)
             call dlnydlnx(log(Mint),log(sigma),Mres,dlnsdlnM)

             !Calculate number of clusters on mass grid points (observed mass)
             do j = 1,Mres
                cmf = 0.d0

                !Integration over true cluster mass taking into account mass scattering (eq. 4.3)
                do k = 1,Mres

                   !Eq. 4.3 for a mass bin of width dM
                   w = MassScatter(Mint(k),Mint(j),Mmin,Mint(j)*dlogM)
                   if (k.eq.1 .or. k.eq.Mres) w = 0.5d0*w

                   !Calculating the cluster mass function, ie. the number of clusters of a given true mass (dn/dm)
                   hmf = clustermassfunc(Mint(k),sigma(k),dlnsdlnM(k),delta_c)

                   !w*hmf*dM is the number of clusters of true mass Mint(k) detected to have mass Mint(j)
                   cmf = cmf + w*hmf*Mint(k)*dlogM
                end do

                !Multiply the cluster mass function with comoving volume, dz, and sky coverage
                pop(i,j) = cmf*CP%fskycl*4.d0*pi*dz*dV
             end do

             !If detection threshold exceeds Mhi (10^16 Msun)
          else

             !Find the lowest redshift where detection threshold exceeds Mhi
             if (z .lt. zhi) then
                zhi = z
             end if

             !No clusters if detection threshold exceeds Mhi
             do j = 1,Mres
                pop(i,j) = 0.d0
             end do
          end if

       else
          do j = 1,Mres
             pop(i,j) = 0.d0
          end do
       end if
       
    end do

    !$OMP END PARAllEl DO

  end subroutine populate

  !Function for computing the mass scatter before binning the clusters, ie. a mass bin with a single point and width dM
  function MassScatter(M,Mobs,Mmin,dM)
    implicit none
    real(dl) :: M,Mmin,MassScatter,lnM_mean,Mobs,dM

    lnM_mean = log(M)-CP%lnM_sigma**2.d0/2.d0

    MassScatter = 0.d0

    if(Mobs.ge.Mmin) MassScatter = exp(-(log(Mobs)-lnM_mean)**2.d0/2.d0/CP%lnM_sigma**2.d0)/(Mobs*CP%lnM_sigma*sqrt(2*pi))*dM

    return
  end function MassScatter

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routines for filling up redshift and mass                 %%
  !%%           with approx. equal numbers of clusters in                 %%
  !%%           each bin. This is the routines called when                %%
  !%%           code is run from CosmoMC on a generated                   %%
  !%%           dataset.                                                  %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !Routine for filling redshift and mass bins (selected based on fiducial model).
  subroutine FillClusterBins(clustersperbin,zb,Mb)
    implicit none
    real(dl) :: clustersperbin(CP%cl_zbins,CP%cl_Mbins),zb(cl_zbins_max+1),Mb(cl_zbins_max,cl_Mbins_max+1)
    integer :: i,j,l

    clustersperbin = 0.d0

    !$OMP PARALLEL DO DEFAUlT(SHARED) SCHEDULE(DYNAMIC) &
    !$OMP & PRIVATE(i,j,l)
    do l = CP%cl_zbins*CP%cl_Mbins-1,0,-1
       i = (l-mod(l,CP%cl_Mbins))/CP%cl_Mbins+1
       j = mod(l,CP%cl_Mbins)+1
       clustersperbin(i,j) = rombint_obj_mod((/Mb(i,j),Mb(i,j+1),0.d0,0.d0,0.d0/),ClustersInZnMRange,zb(i),zb(i+1),tol_clustersperbin,20)
    end do

    !$OMP END PARALLEL DO

    return
  end subroutine FillClusterBins

  function ClustersInZnMRange(obj,z)
     implicit none
     real(dl) :: z,delta_c,Mmin,dV,cmf,Delta_vir,Mbl,Mbh
     real(dl) :: intlo,inthi,s2
     real(dl), dimension(5) :: obj
     real(dl) :: ClustersInZnMRange

     Mbl = obj(1)
     Mbh = obj(2)
     delta_c = 0.d0
     Delta_vir = 0.d0
     s2 = CP%lnM_sigma**2.d0

     !For each redshift value the spherical collapse is solved and delta_c(z) and Delta_vir(z) are returned
     if (z .ne. 0.) call spherical_collapse(z,delta_c,Delta_vir) 

     if (delta_c .gt. 0.d0 .and. Delta_vir .gt. 0.d0 .and. z .ne. 0.) then

        !Computation of the detectiong threshold at z
        call rtsafe(CC_Mmin,z,Mlo,Mhi,acc_Mmin,Mmin,Delta_vir) !calculate Mmin at given redshift !usikkerhed på 1d3 OEB

        !If the detection threshold does not exceed Mhi(10^16 Msun) - integration over Dark Matter mass is preformed
        if ((Mmin.lt.Mhi).and.(Mmin.gt.0.d0))then !linear integration
              !Narrow down integration range to region with contribution to save time
              intlo = log(Mbl)-sqrt(-2.0*s2*log(Mbl/Mhi*sqrt(2.d0*pi*s2)*1.d-6))+s2/2.d0
              inthi = log(Mbh)+sqrt(-2.0*s2*log(Mbh/Mhi*sqrt(2.d0*pi*s2)*1.d-6))+s2/2.d0

              cmf = rombint_obj_mod((/Mbl,Mbh,z,Mmin,delta_c/),ClustersInMassRange,intlo,inthi,tol_cmf,20)
              !Comoving volume element for eq. 3.7 
              dV = comvol(z)

              !Adding to the trapezoidal calculation of the redshift integral
              ClustersInZnMRange = cmf*CP%fskycl*4.d0*pi*dV
        end if
     else
        ClustersInZnMRange = 0.d0
     end if
     return
  end function ClustersInZnMRange

  function ClustersInMassRange(obj,lnM)
     implicit none
     real(dl) :: ClustersInMassRange
     real(dl), dimension(5) :: obj
     real(dl) :: M,lnM,z,Mmin,f,s,dlnsdlnM,w,delta_c,Mbl,Mbh

     Mbl = obj(1)
     Mbh = obj(2)
     z = obj(3)
     Mmin = obj(4)
     delta_c = obj(5)

     M = exp(lnM)
     call GSigma(M,z,s,dlnsdlnM)
     dlnsdlnM = dlnsdlnM*M/s
     w = (erf((log(max(Mmin,Mbh))-lnM+CP%lnM_sigma**2.d0/2.d0)/sqrt(2.d0)/CP%lnM_sigma)&
         -erf((log(max(Mmin,Mbl))-lnM+CP%lnM_sigma**2.d0/2.d0)/sqrt(2.d0)/CP%lnM_sigma))/2.d0
     f = clustermassfunc(M,s,dlnsdlnM,delta_c)
     ClustersInMassRange = w*f*M
     return
  end function ClustersInMassRange

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Sheth-Tormen cluster mass function presented              %%
  !%%           in astro-ph/9901122. Edit to use your own                 %%
  !%%           mass function.                                            %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Fitting cluster mass function, we use Sheth-Tormen astro-ph/9901122
  function clustermassfunc(M,S,dlnsdlnM,delta_c)
    implicit none
    real(dl) :: M,clustermassfunc,dndM_a,dndM_X,dndM_p,rhobar,S,dlnSdlnM,delta_c
    real(dl), parameter :: Mpl2Msun=9.1385d37
    real(dl), parameter :: uniconv=2.9979d8**2*3.0857d22/1.989d30/(8._dl*pi*6.673d-11)!grhom2Msun/Mpc³

    !Fitting parameters for Sheth-Tormen
    dndM_a=0.707
    dndM_X=0.322184
    dndM_p=0.3

    !Present background matter density
    rhobar=(CP%omegac+CP%omegab)*grhom*uniconv

    !Sheth-Tormen - see eq. 3.1
    clustermassfunc=-sqrt(2._dl*dndM_a/pi)*dndM_X*(1._dl+(dndM_a*delta_c**2/S**2)**(-dndM_p))*rhobar/M**2*delta_c/S*dlnSdlnM*exp(-dndM_a*delta_c**2/2._dl/S**2) !for linear-integration

    return
  end function clustermassfunc

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routine for calculating the variance of the               %%
  !%%           matter density field.                                     %%
  !%%                                                                     %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !Calculate the variance of the matter density field and derivative with
  !respect to M
  subroutine GSigma(M,z,S,dSdM)
    implicit none
    real(dl) :: M,R,rhobar,z,transf,h,a,b,S,hubble,kh,k,x,win,dS,dSo,lnk,lnko,dlnk,delta
    real(dl) :: dSdM,ddSdM,ddSdMo,dwin,cosx,sinx,dRdM,x2
    integer :: zhi,zlo,ztemp,ik
    real(dl), parameter :: uniconv=2.9979d8**2*3.0857d22/1.989d30/(8._dl*pi*6.673d-11)!grhom2Msun/Mpc³

    !Present background matter density (Msun/Mpc³)
    rhobar=(CP%omegac+CP%omegab)*grhom*uniconv  

    !Comoving smoothing scale
    R=(3._dl*M/(4._dl*pi*rhobar))**(1._dl/3._dl) 
    dRdM = (3.d0/4.d0/pi/rhobar)**(1.d0/3.d0)*M**(-2.d0/3.d0)/3.d0

    !Redshift indices
    zhi = CS_num_z_tf
    zlo = 1

    !Hubble constant
    hubble = CP%h0/100.d0

    !Find indices of redshift values in table of transfer functions that bracket relevant redshift
    do while (zhi-zlo .gt. 1)
       ztemp=(zhi+zlo)/2

       if(CP%Transfer%redshifts(ztemp) .lt. z) then
          zhi=ztemp
       else
          zlo=ztemp
       endif
    end do

    !Linear interpolation of transfer function
    h=CP%Transfer%redshifts(zhi)-CP%Transfer%redshifts(zlo)

    a=(CP%Transfer%redshifts(zhi)-z)/h
    b=(z-CP%Transfer%redshifts(zlo))/h           

    !Trapezoidal integral over ln(k)
    S = 0.d0
    dSo = 0.d0
    dSdM = 0.d0
    ddSdMo = 0.d0
    do ik=1,num_q_trans
       kh = q_trans(ik)
       if (kh == 0) cycle
       k = kh*hubble
       delta = k**2.d0*(a*LinTransferData(ik,zlo)+b*LinTransferData(ik,zhi))

       x = k*R
       x2 = x*x
       sinx = sin(x)
       cosx = cos(x)
       win = 3.d0*(sinx-x*cosx)/x2/x
       dwin = 3.d0*(sinx/x2-3.d0*(sinx/x2/x2-cosx/x2/x))
       lnk = log(k)

       if(ik == 1) then
          dlnk = 0.5d0
       else
          dlnk = lnk-lnko
       end if

       dS = (win*delta)**2.d0*ScalarPower(k,1)
       S = S+(dS+dSo)*dlnk/2.d0
       dSo = dS
      
       ddSdM = 2*win*dwin*delta**2.d0*k*ScalarPower(k,1)
       dSdM = dSdM+(ddSdM+ddSdMo)*dlnk/2.d0
       ddSdMo = ddSdM

       lnko = lnk
    end do

    S = sqrt(S)
    dSdM = dSdM/2.d0/S*dRdM
    return
  end subroutine GSigma

  !Calculate the variance of the matter density field
  function Get_sigma(M,z)
    implicit none
    real(dl) :: M,R,rhobar,z,transf,h,a,b,Get_sigma,S,hubble,kh,k,x,win,dS,dSo,lnk,lnko,dlnk,delta
    integer :: zhi,zlo,ztemp,ik
    real(dl), parameter :: uniconv=2.9979d8**2*3.0857d22/1.989d30/(8._dl*pi*6.673d-11)!grhom2Msun/Mpc³

    !Present background matter density (Msun/Mpc³)
    rhobar=(CP%omegac+CP%omegab)*grhom*uniconv  

    !Comoving smoothing scale
    R=(3._dl*M/(4._dl*pi*rhobar))**(1._dl/3._dl) 		 

    !Redshift indices
    zhi = CS_num_z_tf
    zlo = 1

    !Hubble constant
    hubble = CP%h0/100.d0

    !Find indices of redshift values in table of transfer functions that bracket relevant redshift
    do while (zhi-zlo .gt. 1)
       ztemp=(zhi+zlo)/2

       if(CP%Transfer%redshifts(ztemp) .lt. z) then
          zhi=ztemp
       else
          zlo=ztemp
       endif
    end do

    !Linear interpolation of transfer function
    h=CP%Transfer%redshifts(zhi)-CP%Transfer%redshifts(zlo)

    a=(CP%Transfer%redshifts(zhi)-z)/h
    b=(z-CP%Transfer%redshifts(zlo))/h           

    !Trapezoidal integral over ln(k)
    S = 0.d0
    do ik=1,num_q_trans
       kh = q_trans(ik)
       if (kh == 0) cycle
       k = kh*hubble
       delta = k**2.d0*(a*LinTransferData(ik,zlo)+b*LinTransferData(ik,zhi))

       x = k*R
       win = 3.d0*(sin(x)-x*cos(x))/x**3.d0
       lnk = log(k)

       if(ik == 1) then
          dlnk = 0.5d0
       else
          dlnk = lnk-lnko
       end if

       dS = (win*delta)**2.d0*ScalarPower(k,1)
       S = S+(dS+dSo)*dlnk/2.d0
       dSo = dS
       lnko = lnk
    end do

    Get_sigma = sqrt(S)
    return
  end function Get_sigma

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routines for determining the detection                    %%
  !%%           threshold of the cluster survey.                          %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !The zero point of this subroutine determines Mmin (detection threshold)
  subroutine CC_Mmin(M,z,out,Delta_vir)
    implicit none
    real(dl) :: z,M,out
    real(dl) :: sigma_noise,kappa_G,R_s,R_vir,rhobar,Delta_vir,a,nbg,rombint
    real(dl), parameter :: Mpl2Msun=9.1385d37
    real(dl), parameter :: uniconv=2.9979d8**2*3.0857d22/1.989d30/(8._dl*pi*6.673d-11)!grhom2Msun/Mpc³
    external :: rombint     

    !Present background matter density
    rhobar=(CP%omegac+CP%omegab)*grhom*uniconv

    !nbg = CP%ngal_cl-rombint(galaxy_basic2,0.d0,z,tol_sigmaCR,20)
    
    !Noise in the shear signal eq. 5.5
    sigma_noise=sqrt(CP%mean_int_ellip_cl**2/(4._dl*pi*CP%theta_G**2*(CP%ngal_cl)))
    !sigma_noise=sqrt(CP%mean_int_ellip_cl**2/(4._dl*pi*CP%theta_G**2*(nbg/3600*206265*206265)))

    !Smoothing of the NFW profile with Gaussian filter eq. 5.2
    a=alpha(z,M,Delta_vir)

    !Minimal shear signal for detection
    kappa_G=CP%signal2noise*sigma_noise

    !Calculate virial radius required for cluster of mass M to give significant signal (comment is not entirely correct since eq. 5.2 (alpha) also depends on R_vir)
    R_s=sqrt(M*a/(SigmaCR(z,CP%ngal_cl)*pi*kappa_G))
    !R_s=sqrt(M*a/(SigmaCR(z,nbg)*pi*kappa_G))
    R_vir=R_s*CP%c_nfw

    !Discrepency between M and calculated R_vir - the root of this gives the detection threshold
    out=4._dl/3._dl*pi*R_vir**3.d0*rhobar*(1.d0+z)**3.d0*(Delta_vir)-M

    return
  end subroutine CC_Mmin

  !Function used in CC_Mmin - average over gaussian aperture - eq. 5.2
  function alpha(z,M,Delta_vir)
    implicit none
    real(dl) :: z,M,alpha,dx,alpha_test,rombint_obj,Delta_vir,alpha_lin
    integer :: i,num_x
    real :: x_G
    external rombint_obj

    x_G = CP%theta_G*pi/10800.d0/theta_s(z,M,Delta_vir)

    !Use the rombin_obj already in CAMB - the integral in log-scale
    alpha = rombint_obj(x_G,alphafunc_ln,log(1.d-10),log(CP%c_nfw),tol_alpha,30)
    alpha=alpha/(log(1._dl+CP%c_nfw)-CP%c_nfw/(1._dl+CP%c_nfw))
    return
  end function alpha

  !Function to be integrated in alpha(z,M) - log-scale
  function alphafunc_ln(x_dummy,lnx)
    implicit none
    real(dl) :: x,alphafunc_ln,x_G,lnx
    real :: x_dummy
    x_G = x_dummy
    x = exp(lnx)
    alphafunc_ln = x*x/x_G**2.d0*exp(-(x/x_G)**2.d0)*f(x)  !f(x) is truncated nfw profile - ArXiv 0310607v2
    return
  end function alphafunc_ln

  !Truncated nfw profile - ArXiv 0310607v2				
  function f(x)
    implicit none
    real(dl) :: x,f
    if(x.lt.1)then
       f=-sqrt(CP%c_nfw**2-x**2)/((1._dl-x**2)*(1._dl+CP%c_nfw))+1._dl/(1._dl-x**2)**(3.d0/2.d0)*acosh((x**2+CP%c_nfw)/(x*(1._dl+CP%c_nfw)))
    elseif(x.eq.1)then
       f=sqrt(CP%c_nfw**2-1._dl)/(3._dl*(1._dl+CP%c_nfw))*(1._dl+1._dl/(1._dl+CP%c_nfw))
    else
       f=-sqrt(CP%c_nfw**2-x**2)/((1._dl-x**2)*(1._dl+CP%c_nfw))-1._dl/(x**2-1._dl)**(3.d0/2.d0)*acos((x**2+CP%c_nfw)/(x*(1._dl+CP%c_nfw)))
    endif
  end function f

  !Function used to calculate alpha - eq. 5.2
  function theta_s(z,M,Delta_vir)
    implicit none
    real(dl) :: z,theta_s,R_s,R_vir,rhobar,M,Delta_vir
    real(dl), parameter :: Mpl2Msun=9.1385d37
    real(dl), parameter :: uniconv=2.9979d8**2*3.0857d22/1.989d30/(8._dl*pi*6.673d-11)!grhom2Msun/Mpc³
    rhobar=(CP%omegac+CP%omegab)*grhom*uniconv
    R_vir=(M/(4._dl/3._dl*pi*rhobar*((1.d0+z)**3.d0)*(Delta_vir)))**(1.d0/3.d0)
    R_s=R_vir/CP%c_nfw
    theta_s=R_s*(1._dl+z)/comdist(z)
    return
  end function theta_s

  !Function integrated by SigamCR(z_l) - eq. 5.3
  function Sigmafunc(z_dummy,z)
    implicit none
    real(dl) :: z,z_l,Sigmafunc
    real :: z_dummy
    z_l = z_dummy
    Sigmafunc=galaxy_basic2(z)*(1._dl-comdist(z_l)/comdist(z))
    return
  end function Sigmafunc

  !Calculates the critical surface density - eq. 5.3
  function SigmaCR(z_l,nbg)
    implicit none
    real(dl) :: z_l,SigmaCR,z_max,rombint_obj,nbg
    real :: z_dummy
    integer :: i,num_sigma
    external rombint_obj

    z_max=2._dl

    z_dummy = z_l

    sigmaCR = rombint_obj(z_dummy,Sigmafunc,z_l,z_max,tol_sigmaCR,20)

    !Inverse critical surface mass density in Mpc^2/Msun
    SigmaCR=SigmaCR*4._dl*pi*6.673d-11*1.989d30/((2.9979d8**2*3.086d22)*(1._dl+z_l))*comdist(z_l)/nbg

    SigmaCR=1._dl/SigmaCR

    return
  end function SigmaCR


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Routines for solving the spherical collapse.              %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !Subroutine for the spherical collapse  
  subroutine spherical_collapse(z,dcf,Dvf)

    implicit none

    double precision :: z,dcf,Dvf,deltami,deltamig,dummy,dcf0,Dvf0
    integer :: counter

    counter = 1
    deltami = -1.d0

    call set_delta(deltamig,z)	!Guess for initial matter over density from 2D polynomial - for faster convergence (works but might need to be updated for w0,wa parametrization)
    do while (deltami .le. 0.d0 .and. counter .le. 5)
       call rtsafe(solve,z,0.1d0*deltamig,10.d0*deltamig,acc_deltami,deltami,dummy)  !Determine initial matter over density corresponding to collapse at z 
       counter = counter + 1
    end do

    Dvf = 0.0d0
    if (deltami .gt. 0.d0) then
       call Get_deltas(deltami,z,dcf,Dvf)  !Calculate delta_c and Delta_vir with the right initial matter over density
    else 
       dcf = 0.d0
    end if

    return
  end subroutine spherical_collapse

  !Differential equations for spherical collapse
  subroutine diffeq(tau,y,dydt,deltami,Xi)
    implicit none
    double precision, dimension(5),intent(in) :: y
    double precision :: tau,om,deltami,Xi,om0,om0_tot,om_tot
    double precision, dimension(5) :: dydt

    om0_tot = CP%omegab+CP%omegac+CP%omegan
    om0 = CP%omegab+CP%omegac
    call omegaM(y(3),om0,om)
    call omegaM(y(3),om0_tot,om_tot)

    dydt(1) = y(2);  !Xprime
    dydt(2) = -2.0d0*sqrt(om0_tot/y(3)/om_tot)/y(3)*y(2)-0.5d0*(om0/(y(3)**3.0d0)*(((Xi/y(1))**3.0d0)*(1.0d0+deltami)-1.0d0))*y(1)  !Xprimeprime       
    dydt(3) = sqrt(om0_tot/y(3)/om_tot);  !aprime/a0
    dydt(4) = -1.0d0/y(3)*y(5)/(CP%H0/1000.d0);    !delta_m lin prime
    dydt(5) = -sqrt(om0_tot/y(3)/om_tot)/y(3)*y(5)-3.0d0/2.0d0*(om0/(y(3)**2.0d0)/om)*(CP%H0/1000.d0)*(om*y(4));   !theta_m lin prime
    return
  end subroutine diffeq

  !Solves the spherical collapse - root of return value newtout determines the right initial matter over density
  subroutine solve(deltami,z,newtout,dummy)
    implicit none
    double precision, dimension(5) :: y
    double precision :: tau,h,om_tot,newtout,z,deltami,Xi,Dvir,om0_tot,om0,ai,dummy

    h = 1.d-10
    om0_tot = CP%omegac+CP%omegab+CP%omegan
    om0 = CP%omegac+CP%omegab
    ai = (3.d0*CP%cl_taui*sqrt(om0_tot)/2.d0)**(2.d0/3.d0)
    call omegaM(ai,om0_tot,om_tot)
    Xi = (2.0d0*G*CP%cl_M/(CP%H0/1000.d0)**2.0d0/om0*(3.0857d0/1.989d5)/(1.d0+deltami))**(1.0d0/3.0d0)
    y(1) = Xi;			!comoving radius of sphere X
    y(2) = (2.0d0/3.0d0/CP%cl_taui*(1.0d0-deltami/3.0d0)-sqrt(om0_tot/ai/om_tot)/ai)*Xi;		!Xprime derivative of X - Xprime = Rprime/R-aprime
    y(3) = ai;		!a/a0
    y(4) = deltami;	!delta_m_linear
    y(5) = -ai*deltami*sqrt(om0_tot/ai/om_tot)/ai*CP%H0/1000.d0;

    call odeint(y,5,CP%cl_taui,10.d0,eps_ode,h,1.d-10,diffeq,rkqs,deltami,Xi,Dvir,z,.false.)

    newtout = 1.d0/y(3)-1.d0-z 
    return
  end subroutine solve

  !Same as solve, but returns deltas c_s -> inf
  subroutine Get_deltas(deltami,z,dcrit,Dvir)
    implicit none
    double precision, dimension(5) :: y
    double precision :: tau,h,om_tot,newtout,dcrit,Dvir,deltami,Xi,om0,om0_tot,ai,z,xend
    logical :: doVir
    xend = 10.d0
    h = 1.d-10
    om0_tot = CP%omegac+CP%omegab+CP%omegan
    om0 = CP%omegac+CP%omegab
    ai = (3.d0*CP%cl_taui*sqrt(om0_tot)/2.d0)**(2.d0/3.d0)
    call omegaM(ai,om0_tot,om_tot)
    Xi = (2.0d0*G*CP%cl_M/(CP%H0/1000.d0)**2.0d0/om0*(3.0857d0/1.989d5)/(1.d0+deltami))**(1.0d0/3.0d0)

    y(1) = Xi;			!comoving radius of sphere X
    y(2) = (2.0d0/3.0d0/CP%cl_taui*(1.0d0-deltami/3.0d0)-sqrt(om0_tot/ai/om_tot)/ai)*Xi;		!Xprime derivative of X - Xprime = Rprime/R-aprime
    y(3) = ai;		!a/a0
    y(4) = deltami;	!delta_m_linear
    y(5) = -ai*deltami*sqrt(om0_tot/ai/om_tot)/ai*CP%H0/1000.d0;

    call odeint(y,5,CP%cl_taui,xend,eps_ode,h,1.d-10,diffeq,rkqs,deltami,Xi,Dvir,z,.false.)
    xend = Dvir
    dcrit = y(4)
    Dvir = 0.d0
    call odeint(y,5,xend,CP%cl_taui,eps_ode,h,1.d-10,diffeq,rkqs,deltami,Xi,Dvir,z,.true.)
    return
  end subroutine Get_deltas

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Ordinary differential equation solver                     %%
  !%%           from Numerical recipes with adaptive                      %%
  !%%           step size and Cash-Karp Runge Kutta                       %%
  !%%           stepping. Modified to suit need for                       %%
  !%%           calculation of virialisation.                             %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! This routine is adapted from numerical recipes and modified for virialisation calculation
  subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,derivs,rkqs,deltami,Xi,Dvir,z,doDvir)
    implicit none
    integer :: nvar,MAXSTP
    logical :: doDvir
    double precision :: eps,h1,hmin,x1,x2,ystart(nvar),TINY,deltami,Xi,Dvir,z
    external derivs, rkqs
    parameter (MAXSTP=10000,TINY=1.d-30)
    integer :: i,kmax,kount,nstp
    double precision :: dxsav,h,hdid,hnext,x,xsav,appa,apa,xppx,xpx,vir,om,eos,xold,om0,om0_tot,om_tot
    double precision :: dydx(nvar),y(nvar),yscal(nvar),yold(nvar)

    kmax = 1000
    dxsav = 1.d-8

    x = x1
    h = sign(h1,x2-x1)
    kount = 0
    do i = 1,nvar
       y(i) = ystart(i)
    enddo

    if (kmax .gt. 0) xsav = x-2.d0*dxsav

    do nstp = 1,MAXSTP

       call derivs(x,y,dydx,deltami,Xi)

       do i = 1,nvar
          yscal(i) = abs(y(i)) + abs(h*dydx(i)) + TINY
       enddo
       if (kmax .gt. 0) then
          if (abs(x-xsav) .gt. abs(dxsav)) then
             if(kount .lt. kmax-1) then
                kount = kount + 1
                xsav = x
             endif
          endif
       endif

       if ((x+h-x2)*(x+h-x1) .gt. 0.d0) h = x2-x

       yold = y
       xold = x

       call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,deltami,Xi)

       !The following takes care of Delta_vir requiring that the second derivative of the moment of inertia with respect to cosmic time is zero
       if ((doDvir .eqv. .true.)) then 
          om0_tot = CP%omegac+CP%omegab+CP%omegan
          om0 = CP%omegac+CP%omegab
          call omegaM(y(3),om0,om)
          call omegaM(y(3),om0_tot,om_tot)
          eos = CP%w0+(1.d0-y(3))*CP%wa

          appa = -1.d0/2.d0*(om0*y(3)**(-3.d0)+(1.d0-om0)*exp(-3.d0*CP%wa*(1.d0-y(3)))*((1.d0+3.d0*eos)*y(3)**(-3.d0*(1.d0+CP%w0+CP%wa))))
          apa = sqrt(om0_tot/y(3)/om_tot)/y(3)
          xpx = y(2)/y(1)

          xppx = (-2.d0*sqrt(om0_tot/y(3)/om_tot)/y(3)*y(2)-0.5d0*om0/y(3)**3.d0/om*(om*((Xi/y(1))**3.d0*(1.d0+deltami)-1.d0))*y(1))/y(1)
          vir = xppx+xpx**2.d0+apa**2.d0+4.d0*apa*xpx+appa

          if ((Dvir .eq. 0.d0) .and. (abs(vir) .le. 1.d-4)) then
             Dvir = (1.d0+deltami)*(Xi/y(1))**3.d0/(1.d0+z)**3.d0/y(3)**3.d0
             return
          endif
          if ((Dvir .eq. 0.d0) .and. (vir .lt. -1.d-4)) then   !redo step
             hnext = hnext/10
             y = yold
             x = xold
          endif

       endif
       if (doDvir .eqv. .false.) then
          if ((y(1) .lt. 1.d-2) .or. (y(3) .gt. 2.d0)) then
             do i = 1,nvar
                ystart(i) = y(i)
             enddo
             if (kmax .ne. 0) then
                kount = kount + 1
             endif
             Dvir = x

             return
          endif
       endif

       if ((x-x2)*(x2-x1) .ge. 0.d0) then
          do i = 1,nvar
             ystart(i) = y(i)
          enddo
          if (kmax .ne. 0) then
             kount = kount + 1
          endif

          return
       endif

       h = hnext
    enddo

    return
  end subroutine odeint

  ! This routine is adapted from numerical recipes.
  subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs,deltami,Xi)
    implicit none
    integer :: n
    double precision :: eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),deltami,Xi
    external derivs

    integer :: i
    double precision :: errmax,h,htemp,yerr(n),ytemp(n),SAFETY,PGROW,PSHRNK,ERRCON
    parameter (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)

    h = htry

1   call rkck(y,dydx,n,x,h,ytemp,yerr,derivs,deltami,Xi)

    errmax = 0.d0

    do i = 1,n
       errmax = max(errmax,abs(yerr(i)/yscal(i)))
    enddo

    errmax = errmax/eps

    if (errmax .gt. 1.d0) then
       htemp = SAFETY*h*(errmax**PSHRNK)
       h = sign(max(abs(htemp),0.1d0*abs(h)),h)
       goto 1
    else
       if(errmax .gt. ERRCON) then
          hnext = SAFETY*h*(errmax**PGROW)
       else
          hnext = 5.d0*h
       endif
       hdid = h
       x = x+h
       do i = 1,n
          y(i) = ytemp(i)
       enddo
       return
    endif

  end subroutine rkqs

  ! This routine is adapted from numerical recipes.
  subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs,deltami,Xi)
    implicit none
    integer :: n
    double precision :: h,x,dydx(n),y(n),yerr(n),yout(n),deltami,Xi
    external derivs

    integer i
    double precision, DIMENSION(n) :: ak2,ak3,ak4,ak5,ak6,ytemp
    double precision, PARAMETER :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,&
         A6=0.875d0,B21=0.2d0,B31=3.0d0/40.0d0,B32=9.0d0/40.0d0,&
         B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.0d0/54.0d0,&
         B52=2.5d0,B53=-70.0d0/27.0d0,B54=35.0d0/27.0d0,&
         B61=1631.0d0/55296.0d0,B62=175.0d0/512.0d0,&
         B63=575.0d0/13824.0d0,B64=44275.0d0/110592.0d0,&
         B65=253.0d0/4096.0d0,C1=37.0d0/378.0d0,&
         C3=250.0d0/621.0d0,C4=125.0d0/594.0d0,&
         C6=512.0d0/1771.0d0,DC1=C1-2825.0d0/27648.0d0,&
         DC3=C3-18575.0d0/48384.0d0,DC4=C4-13525.0d0/55296.0d0,&
         DC5=-277.0d0/14336.0d0,DC6=C6-0.25d0


    do i = 1,n
       ytemp(i) = y(i) + B21*h*dydx(i)
    enddo

    call derivs(x+A2*h,ytemp,ak2,deltami,Xi)

    do i = 1,n
       ytemp(i) = y(i) + h*(B31*dydx(i)+B32*ak2(i))
    enddo

    call derivs(x+A3*h,ytemp,ak3,deltami,Xi)

    do i = 1,n
       ytemp(i) = y(i) + h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    enddo

    call derivs(x+A4*h,ytemp,ak4,deltami,Xi)

    do i = 1,n
       ytemp(i) = y(i) + h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
    enddo

    call derivs(x+A5*h,ytemp,ak5,deltami,Xi)

    do i = 1,n
       ytemp(i) = y(i) + h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    enddo

    call derivs(x+A6*h,ytemp,ak6,deltami,Xi)

    do i = 1,n
       yout(i) = y(i) + h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
    enddo

    do i = 1,n
       yerr(i) = h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
    enddo

    return
  end subroutine rkck

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Root finding algorithm based on bisection                 %%
  !%%           and Newton-Raphson from Numerical Recipes.                %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! This routine is adapted from numerical recipes - modified with program
  ! specific exits
  subroutine rtsafe(funcd,z,x1,x2,xacc,xroot,Delta_vir)
    implicit none
    integer :: MAXIT
    double precision :: xroot,x1,x2,xacc,Delta_vir
    external funcd
    parameter (MAXIT=100)

    integer :: j
    double precision :: df,dx,dxold,f,fh,fl,temp,xh,xl,ftemp,z,xt

    call funcd(x1,z,fl,Delta_vir)
    call funcd(x2,z,fh,Delta_vir)

    if (x1 .eq. 1.d0) then   
       if (((fl .lt. 0.d0) .and. (fh .lt. 0.d0)) .or. ((fl .gt. 0.d0) .and. (fh .gt. 0.d0))) then
          xroot = -10.d0
       else
          xroot = -1.d0
       end if
       return
    end if
    
32  if (((fl .lt. 0.d0) .and. (fh .lt. 0.d0)) .or. ((fl .gt. 0.d0) .and. (fh .gt. 0.d0))) then
       xroot = -10.d0
       return
    else if (x2 .lt. 1.d0 .and. fl .lt. -0.5d0) then
       xt = 0.5d0*(x1+x2)
31     call funcd(xt,z,fl,Delta_vir)
       if (((fl .lt. 0.d0) .and. (fh .lt. 0.d0)) .or. ((fl .gt. 0.d0) .and. (fh .gt. 0.d0))) then
          xt = 0.5d0*(x1+xt)
          goto 31
       endif
       x1 = xt
       goto 32
    else if (fl .eq. 0.d0) then
       xroot = x1
       return
    else if (fh .eq. 0.d0) then
       xroot = x2
       return
    else if (fl .lt. 0.d0) then
       xl = x1
       xh = x2
    else
       xh = x1
       xl = x2
    endif

    xroot = .5d0*(x1+x2)
    dxold = abs(x2-x1)
    dx = dxold

    call funcd(xroot,z,f,Delta_vir)
    call funcd(xroot*(1.d0+1.d-4),z,ftemp,Delta_vir)

    df = (ftemp-f)/xroot/1.d-4

    do j = 1,MAXIT
       if (((xroot-xh)*df-f)*((xroot-xl)*df-f) .gt. 0.d0 .or. abs(2.d0*f) .gt. abs(dxold*df)) then
          dxold = dx
          dx = .5d0*(xh-xl)
          xroot = xl+dx
          if (xl .eq. xroot) return
       else
          dxold = dx
          dx = f/df
          temp = xroot
          xroot = xroot-dx
          if (temp .eq. xroot) return
       endif

       if (abs(dx) .lt. xacc) return

       call funcd(xroot,z,f,Delta_vir)
       call funcd(xroot*(1.d0+1.d-4),z,ftemp,Delta_vir)

       df = (ftemp-f)/xroot/1.d-4

       if (f .lt. 0) then
          xl = xroot
       else
          xh = xroot
       endif
    enddo

    return
  end subroutine rtsafe

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Smaller utility functions.                                %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !Ln representation of arccosh
  function acosh(x)
    implicit none
    real(dl) :: x,acosh
    acosh=log(x+sqrt(x**2-1._dl))
  end function acosh

  !Sets initial guess at matter over density for collapse at redshift ~ z
  subroutine set_delta(deltamig,z)
    implicit none
    double precision p00,p10,p01,p20,p11,p02,deltamig,z
    p00 =     0.00052d0
    p10 =   0.0003562d0
    p01 =    0.000154d0
    p20 =   0.0001309d0
    p11 =  -1.346d-05
    p02 =   1.499d-05
    deltamig = p00 + p10*CP%w0 + p01*z + p20*CP%w0**2.d0 + p11*CP%w0*z + p02*z**2.d0;

    return 
  end subroutine set_delta

  !Omega matter as a function of a
  subroutine omegaM(a,om0,om)
    implicit none
    double precision :: om,a,om0_tot,om0

    om0_tot = CP%omegac+CP%omegab+CP%omegan

    om =  om0/(om0_tot+(1.d0-om0_tot)*a**(-3.0d0*(CP%w0+CP%wa))*exp(-3.d0*CP%wa*(1.d0-a)))
    return
  end subroutine omegaM

  !Calculates comoving distance from dtauda in camb
  function comdist(z)
    implicit none
    real(dl) :: z,comdist,a,rombint,dtauda
    external rombint,dtauda
    a=1._dl/(1._dl+z)
    comdist=rombint(dtauda,a,1._dl,tol)
    return
  end function comdist

  !Calculates comoving volume element for flat cosmology
  function comvol(z)
    implicit none
    real(dl) :: z,a,comvol,dtauda
    external dtauda
    a = 1.d0/(1.d0+z)
    comvol = a**2.d0*comdist(z)**2*dtauda(a)
    return  
  end function comvol

  !Gamma function - from numerical recipes
  function gammln(xx)
    real(dl) :: gammln,xx
    integer :: j
    real(dl) :: ser,tmp,x,y
    real(dl) :: stp = 2.5066282746310005d0	
    real(dl), dimension(6) :: cof = (/76.18009172947146d0,&
         -86.50532032941677d0,24.01409824083091d0,&
         -1.231739572450155d0,0.1208650973866179d-2,&
         -0.5395239384953d-5/)

    x = xx
    y = x
    tmp = x+5.5d0
    tmp = (x+0.5d0)*log(tmp)-tmp
    ser = 1.000000000190015d0
    do j = 1,6
       y = y+1.d0
       ser = ser+cof(j)/y
    enddo

    gammln = tmp+log(stp*ser/x)
    return
  end function gammln

  !Function for computing dln(y)/dln(x) using finite difference
  subroutine dlnydlnx(x,y,n,dy)
    implicit none
    integer :: n,i
    real(dl) :: h,x(n),y(n),dy(n)

    h = x(2)-x(1)
    dy(1) = (-49.0/20.0*y(1)+6.0*y(2)-15.0/2.0*y(3)+20.0/3.0*y(4)-15.0/4.0*y(5)+6.0/5.0*y(6)-1.0/6.0*y(7))/h
    dy(2) = (-49.0/20.0*y(2)+6.0*y(3)-15.0/2.0*y(4)+20.0/3.0*y(5)-15.0/4.0*y(6)+6.0/5.0*y(7)-1.0/6.0*y(8))/h
    dy(3) = (-49.0/20.0*y(3)+6.0*y(4)-15.0/2.0*y(5)+20.0/3.0*y(6)-15.0/4.0*y(7)+6.0/5.0*y(8)-1.0/6.0*y(9))/h
    dy(4) = (-49.0/20.0*y(4)+6.0*y(5)-15.0/2.0*y(6)+20.0/3.0*y(7)-15.0/4.0*y(8)+6.0/5.0*y(9)-1.0/6.0*y(10))/h
    do i = 5,n-4
       dy(i) = (y(i-4)/280.0-4.0*y(i-3)/105.0+y(i-2)/5.0-4.0*y(i-1)/5.0+4.0*y(i+1)/5.0-y(i+2)/5.0+4.0*y(i+3)/105.0-y(i+4)/280.0)/h
    end do
    dy(n-3) = (49.0/20.0*y(n-3)-6.0*y(n-4)+15.0/2.0*y(n-5)-20.0/3.0*y(n-6)+15.0/4.0*y(n-7)-6.0/5.0*y(n-8)+1.0/6.0*y(n-9))/h
    dy(n-2) = (49.0/20.0*y(n-2)-6.0*y(n-3)+15.0/2.0*y(n-4)-20.0/3.0*y(n-5)+15.0/4.0*y(n-6)-6.0/5.0*y(n-7)+1.0/6.0*y(n-8))/h
    dy(n-1) = (49.0/20.0*y(n-1)-6.0*y(n-2)+15.0/2.0*y(n-3)-20.0/3.0*y(n-4)+15.0/4.0*y(n-5)-6.0/5.0*y(n-6)+1.0/6.0*y(n-7))/h
    dy(n) = (49.0/20.0*y(n)-6.0*y(n-1)+15.0/2.0*y(n-2)-20.0/3.0*y(n-3)+15.0/4.0*y(n-4)-6.0/5.0*y(n-5)+1.0/6.0*y(n-6))/h

    !dy = dy*x/y
  end subroutine dlnydlnx

  ! Basic galaxy distribution (this is ultimately the same as galaxy_basic, just normalized)
  function galaxy_basic2(z)
    implicit none
    real(dl) :: z,galaxy_basic2

    galaxy_basic2 = CP%ngal_cl*CP%nz_beta_cl/(CP%nz_z0_cl*exp(gammln((1._dl+CP%nz_alpha_cl)/CP%nz_beta_cl)))*(z/CP%nz_z0_cl)**CP%nz_alpha_cl*exp(-(z/CP%nz_z0_cl)**CP%nz_beta_cl)

    return
  end function galaxy_basic2

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%									  %%
  !%%		Function for outputting CosmoMC data files                %%
  !%%								          %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine OutputSimDataCl(P,OutClusters)
    implicit none

    type (CAMBParams) :: P
    type (OutputClusters) :: OutClusters

    character(LEN=100) :: setoutname,clusteroutname,massbinoutname,datname,massname
    integer :: i

    datname = trim(setname)//"_cluster.dat"
    massname = trim(setname)//"_cluster.massbins"
    setoutname = trim(datarep)//trim(setname)//"_cluster.dataset"
    clusteroutname = trim(datarep)//trim(datname)
    massbinoutname = trim(datarep)//trim(massname)

    call CreateTxtFile(setoutname,37)
    write (37,'(A,A)')  "name = ",trim(setname)
    write (37,'(A,I5)') "sim_ntomobin_cl = ",P%cl_zbins
    write (37,'(A,I5)') "sim_nmassbin_cl = ",P%cl_Mbins
    write (37,'(A,100F15.8)') "zlow_vector_cl =",P%cl_RedshiftBins(1:P%cl_zbins)
    write (37,'(A,100F15.8)') "zhigh_vector_cl =",P%cl_RedshiftBins(2:P%cl_zbins+1)
    write (37,'(A,F18.8)') "efficiency = ",P%efficiency
    write (37,'(A,F18.5)') "completeness = ",P%completeness 
    write (37,'(A,F18.8)') "ngal_cl = ",P%ngal_cl
    write (37,'(A,F18.5)') "nz_z0_cl = ",P%nz_z0_cl
    write (37,'(A,F18.5)') "nz_alpha_cl = ",P%nz_alpha_cl        
    write (37,'(A,F18.5)') "nz_beta_cl = ",p%nz_beta_cl
    write (37,'(A,F18.5)') "c_nfw = ",P%c_nfw
    write (37,'(A,F12.5)') "fskycl = ",P%fskycl
    write (37,'(A,F12.5)') "theta_G = ",P%theta_G
    write (37,'(A,F12.5)') "mean_int_ellip_cl = ",P%mean_int_ellip_cl
    write (37,'(A,F12.5)') "lnM_sigma = ",P%lnM_sigma
    write (37,'(A,F12.5)') "signal2noise = ",P%signal2noise
    write (37,'(A,A)') "datafile = ",trim(datname)
    write (37,'(A,A)') "massbinfile = ",trim(massname)
    close(37)

    call CreateTxtFile(clusteroutname,38)
    do i=1,P%cl_zbins
            write(38,'(1000F18.8)') OutClusters%clustercount(i,1:P%cl_Mbins)
    end do
    close(38)

    call CreateTxtFile(massbinoutname,39)
    do i=1,P%cl_zbins
            write(39,'(1000E15.8)') P%cl_Massbins(i,1:P%cl_Mbins+1)
    end do
    close(39)
 
  end subroutine OutputSimDataCl


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function rombint_obj_mod(obj,f,a,b,tol, maxit)
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, intent(in), optional :: maxit
        integer :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real(dl), dimension(5) :: obj !dummy
        real(dl) f
        external f
        real(dl) :: rombint_obj_mod
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!
	if (present(maxit)) then
            MaxIter = maxit
        end if
        h=0.5d0*(b-a)
        gmax=h*(f(obj,a)+f(obj,b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(obj,a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint_obj_mod=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint_obj_mod,error, tol
        end if
        
        end function rombint_obj_mod

end module Clusters

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
