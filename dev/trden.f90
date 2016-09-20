Subroutine gengrid(origin,dvx,dvy,dvz,nx,ny,nz,grid)

    IMPLICIT None
! Input Variables
    INTEGER, intent(in) :: nx,ny,nz
	REAL(8), intent(in) :: dvx,dvy,dvz
	REAL(8), DIMENSION(3), intent(in) :: origin
! Internal and output Variables
	INTEGER :: ix,iy,iz,nthrd,OMP_GET_NUM_THREADS
    REAL(8), DIMENSION(nx,ny,nz,3), intent(out) :: grid

do iz=1,nz
	do iy=1,ny
		do ix=1,nx
	        grid(ix,iy,iz,1)=origin(1)+(ix-1)*dvx
		    grid(ix,iy,iz,2)=origin(2)+(iy-1)*dvy
			grid(ix,iy,iz,3)=origin(3)+(iz-1)*dvz
        end do
    end do
end do

end Subroutine gengrid


Subroutine diptrde(trden,grid,dv,n,mu)

    IMPLICIT None
! Input Variables
    INTEGER, intent(in) :: n
	REAL(8), intent(in) :: dv
	REAL(8), DIMENSION(n), intent(in) :: trden
	REAL(8), DIMENSION(n,3), intent(in) :: grid
! Internal and output Variables
	INTEGER :: i,OMP_GET_NUM_THREADS
    REAL(8) :: mux,muy,muz
    REAL(8), DIMENSION(3), intent(out) :: mu

mux=0.d0
muy=0.d0
muz=0.d0

do i=1,n
    mux=mux+grid(i,1)*trden(i)*dv
    muy=muy+grid(i,2)*trden(i)*dv
    muz=muz+grid(i,3)*trden(i)*dv
end do

mu(1)=mux
mu(2)=muy
mu(3)=muz

end Subroutine diptrde


Subroutine couptrde(trdena,grida,dva,trdend,gridd,dvd,thresh,na,nd,couptd)

    IMPLICIT None
! Input Variables
    INTEGER, intent(in) :: na,nd
	REAL(8), intent(in) :: dva,dvd,thresh
	REAL(8), DIMENSION(na), intent(in) :: trdena
	REAL(8), DIMENSION(na,3), intent(in) :: grida
	REAL(8), DIMENSION(nd), intent(in) :: trdend
	REAL(8), DIMENSION(nd,3), intent(in) :: gridd
! Internal and output Variables
    REAL(8), DIMENSION(3) :: r12
    REAL(8) :: autown,r,couptd_local,width
	INTEGER :: i,j,OMP_GET_NUM_THREADS
    REAL(8), intent(out) :: couptd

! Description of the Variables:
! r12 is the distance between two points of the grid
! couptd is the coupling

couptd = 0.
autown = 2.194746d5

width = 1. / dsqrt (2*(dva**(2.0/3.0) + dvd**(2.0/3.0)) )
!$omp parallel private(i,j,r12,r,couptd_local)
couptd_local = 0.
!$omp do
do i=1,na
    if (abs(trdena(i)) < thresh) cycle
    do j=1,nd

        if (abs(trdend(j)) < thresh) cycle
        r12=(grida(i,:)-gridd(j,:))
        r=dsqrt(r12(1)**2.0+r12(2)**2.0+r12(3)**2.0)

        couptd_local=couptd_local+trdena(i)*trdend(j)/r * erf(width*r)

    end do
end do
!$omp end do
!$omp critical
couptd=couptd+couptd_local
!$omp end critical
!$omp end parallel
couptd=couptd*dva*dvd*autown

end Subroutine couptrde
