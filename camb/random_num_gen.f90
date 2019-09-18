!	module Precision
! This is from CAMB's subroutines.f90

! 	integer, parameter :: dl = KIND(1.d0)
! 	integer, parameter :: sp = KIND(1.0)
!
! 	!real(dl), parameter :: pi = 3.1415926535897932384626433832795_dl, twopi=2*pi, fourpi=4*pi
! 	!real(dl), parameter :: sqrt6=2.4494897427831780981972840747059_dl

! 	end module Precision

! This module was originally written by Y. Wong for the FutureData package [astro-ph/0606227]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%									%%
!%%	Multivariate random number generator				%%
!%%									%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	module random_numbers
	use Precision 
	private

	integer, parameter :: nmax=10	! Maximum number of correlated variables

	public get_random_cls,nmax

	contains	

	subroutine get_random_cls(nlo,nhi,lmin,lmax,idum,cov_ext,cls_ext)
! This calculates the mock Cls.
	implicit none
	integer, intent(in) :: nlo,nhi,lmin,lmax
	integer, intent(inout) :: idum
	real(dl), intent(in) :: cov_ext(lmin:lmax,nlo:nhi,nlo:nhi)
	real(dl), intent(out) :: cls_ext(lmin:lmax,nlo:nhi,nlo:nhi)

	integer :: i,j,l,m,realn
	real(dl), dimension(:,:,:), allocatable :: alm_matrix
	real(dl), dimension(:,:,:), allocatable :: cov_internal
	real(dl), dimension(:,:,:), allocatable :: cls_internal

	realn=nhi-nlo+1
	
	allocate(cov_internal(lmin:lmax,1:realn,1:realn))
	allocate(cls_internal(lmin:lmax,1:realn,1:realn))
	allocate(alm_matrix(lmin:lmax,-lmax:lmax,1:realn))

	
	! pass to internal variable
	cov_internal(lmin:lmax,1:realn,1:realn)=cov_ext(lmin:lmax,nlo:nhi,nlo:nhi)
	
	
	call get_alm_matrix(realn,lmin,lmax,cov_internal,idum,alm_matrix)

		
	do l=lmin,lmax
		do i=1,realn
			do j=1,i
				cls_internal(l,i,j)=0.0d0
				
				do m=-l,l
					cls_internal(l,i,j)=cls_internal(l,i,j)+alm_matrix(l,m,i)*alm_matrix(l,m,j)
				end do	

				cls_internal(l,i,j)=1.0d0/real(2*l+1,dl)*cls_internal(l,i,j)
			end do
		end do

		do i=1,realn
			do j=i+1,realn
				cls_internal(l,i,j)=cls_internal(l,j,i)
			end do
		end do
	end do
	
	! pass to external variable
	cls_ext(lmin:lmax,nlo:nhi,nlo:nhi)=cls_internal(lmin:lmax,1:realn,1:realn)

	deallocate(alm_matrix)
	deallocate(cov_internal)
	deallocate(cls_internal)

	end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine get_alm_matrix(n,lmin,lmax,cov,idum,alm_matrix)
! Generate correlated random alm values.
! Input covariance matrix cov(lmin:lmax,1:n,1:n), where lmin/lmax is the
! minumum/maximum multipole one wishes to use for the correlations, n is the number
! of correlated variables (e.g., number of tomoraphy bins).
! Results are returned in alm_matrix(lmin:lmax,-lmax:lmax,1:n)

	implicit none
	integer, intent(in) :: n
	integer, intent(in) :: lmin,lmax
	integer, intent(inout) :: idum
	real(dl), intent(in) :: cov(lmin:lmax,n,n)
	real(dl), intent(out) :: alm_matrix(lmin:lmax,-lmax:lmax,n)

	integer :: i,j,m
	real(dl) :: normdev(nmax)

	do i=lmin,lmax
           do j=1,n
              !if (i > 1600) print*,i, cov(i,1:n,j)
           enddo
		do m=-i,i
			do j=1,n
				normdev(j)=gasdev(idum)
				! Generate Gaussian distributed random numbers
			end do

			call MultiRanNormDev(n,cov(i,:,:),normdev(1:n),alm_matrix(i,m,:))
			! This routine turns the random numbers into correlated
			! alm values according to the covariance matrix.
		end do
	end do


	end subroutine get_alm_matrix



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine MultiRanNormDev(n,cov_l,normdev,alm)
! Multivariate Random Normal Deviate generator!!
! Takes the covariance matrix cov_l(1:n,1:n), where n stands for the
! number of correlate variables (e.g., the number of weak lensing tomography
! bins) at each multipole and a 1-D array of Gaussian distributed random
! numbers normdev(1:n) as input.
! Outputs the correlated random multipole moments alm(1:n) (e.g., alm(i)
! could be the convergence of the ith tomography bin).

	implicit none
	integer, intent(in) :: n
	real(dl), intent(in) :: cov_l(n,n)
	real(dl), intent(in) :: normdev(n)
	real(dl), intent(out) :: alm(n)

	integer :: i,j
	real(dl) :: lmatrix(nmax,nmax)

	call choldc(n,cov_l,lmatrix(1:n,1:n))
	! Use Cholesky decomposition to calculate triangular matrix lmatrix


	! This here is the dot product lmatrix.normdev
	do i=1,n
		alm(i)=0.0d0
		do j=1,n
			alm(i)=alm(i)+lmatrix(i,j)*normdev(j)
		end do
	end do

	end subroutine MultiRanNormDev


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine choldc(n,amatrix,lmatrix)
! Cholesky decomposition. The resulting triangular matrix is returned in
! lmatrix(1:n,1:n) in the bottom triangle.
! Input matrix amatrix(1:n,1:n) is unchanged.

	implicit none
	integer, intent(in) :: n
	real(dl), intent(in) :: amatrix(n,n)
	real(dl), intent(out) :: lmatrix(n,n)

	integer :: i,j,k
	real(dl) :: summ
	real(dl) :: p(nmax)

! Passing amatrix[1..n][1..n] to lmatrix[1..n][1..n].
! amatrix[1..n][1..n] will not be used again in this subroutine
! after this operation.
	do i=1,n
		do j=1,n
			lmatrix(i,j)=amatrix(i,j)
		end do
	end do


! Main part of the Cholesky decomposition.
	do i=1,n
		do j=1,n
			summ=lmatrix(i,j)
			do k=i-1,1,-1
				summ=summ-lmatrix(i,k)*lmatrix(j,k)
			end do

			if(i.eq.j) then
                           !print*,lmatrix(i,j)
                           !print*,summ
                           if(summ.le.0.0d0) then
                              print*, i,j, summ, lmatrix(3,3)
                              print*,lmatrix(3,1)**2, lmatrix(3,2)**2
                              pause 'choldc failed'
                           endif
                           p(i)=sqrt(summ)
          		else
                           lmatrix(j,i)=summ/p(i)
			endif
		end do
	end do



! Fill in the diagonal elements.
	do i=1,n
		lmatrix(i,i)=p(i)
	end do


! Fill in the upper triangular matrix with zeros.
	do i=1,n-1
		do j=i+1,n
			lmatrix(i,j)=0.0d0
		end do
	end do

	end subroutine choldc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%	The following routines are random number generators from		%%
!%5	Numerical Recipes.  Replace them if you have something better.		%%
!%%										%%
!%%	1. gasdev(idum) returns a Gaussian distributed random number.  It	%%
!%%	uses ran1(idum) to generate a random number between 0.0 and 1.0.	%%
!%%										%%
!%%	2. ran1(idum) retruns a uniform random deviate between 0.0 and 1.0.	%%
!%%										%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function gasdev(idum)
	!use Precision
	!use random_numbers, ONLY: ran1
! Returns a normally distributed deviate with zero mean and unit variance.
! Initialise with some negative integer idum; do not change idum
! between successive calls.
! From Numerical Recipes.

	implicit none
      	integer, intent(inout) :: idum

	real(dl) :: gasdev,fac,gset,rsq,v1,v2
	integer iset

	SAVE iset,gset
	DATA iset/0/

	if (idum.lt.0) iset=0

	if (iset.eq.0) then
1       	v1=2.0d0*ran1(idum)-1.0d0
        	v2=2.0d0*ran1(idum)-1.0d0
        	rsq=v1*v1+v2*v2

		if(rsq.ge.1.0d0.or.rsq.eq.0.0d0) goto 1

		fac=sqrt(-2.0d0*log(rsq)/rsq)
        	gset=v1*fac
        	gasdev=v2*fac
        	iset=1
      	else
        	gasdev=gset
        	iset=0
      	endif
      	return
      	end function gasdev


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function ran1(idum)
	!use Precision
! Retruns a uniform random deviate between 0.0 and 1.0.
! Initialise with some negative integer idum; do not change idum
! between successive calls.
! From Numerical Recipes

	implicit none
	integer, intent(inout) :: idum

	real(dl) :: ran1

	integer, parameter :: IA=16807
	integer, parameter :: IM=2147483647
	integer, parameter :: IQ=127773
	integer, parameter :: IR=2836
     	integer, parameter :: NTAB=32
	integer, parameter :: NDIV=1+(IM-1)/NTAB

	real(dl), parameter :: EPS=1.2d-7
	real(dl), parameter :: RNMX=1.0d0-EPS
	real(dl), parameter :: AM=1.0d0/IM

	integer :: j,k,iv(NTAB),iy

	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/

	if (idum.le.0.0d0.or.iy.eq.0) then
        	idum=max(-idum,1)

		do 11 j=NTAB+8,1,-1
          		k=idum/IQ
          		idum=IA*(idum-k*IQ)-IR*k

			if (idum.lt.0) idum=idum+IM
          		if (j.le.NTAB) iv(j)=idum
11      	continue
        	iy=iv(1)
      	endif

      	k=idum/IQ
      	idum=IA*(idum-k*IQ)-IR*k

	if (idum.lt.0) idum=idum+IM

      	j=1+iy/NDIV
      	iy=iv(j)
      	iv(j)=idum
      	ran1=min(AM*iy,RNMX)

   	return
	end function ran1

	end module random_numbers

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

