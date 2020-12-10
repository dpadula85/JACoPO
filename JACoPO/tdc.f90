Subroutine coultdc(trdena,grida,dva,trdend,gridd,dvd,thresh,na,nd,couptd)

IMPLICIT None
! Input Variables
INTEGER, intent(in) :: na,nd
REAL(8), intent(in) :: dva,dvd,thresh
REAL(8), DIMENSION(na), intent(in) :: trdena
REAL(8), DIMENSION(nd), intent(in) :: trdend
REAL(8), DIMENSION(na,3), intent(in) :: grida
REAL(8), DIMENSION(nd,3), intent(in) :: gridd
! Internal and output Variables
REAL(8), DIMENSION(3) :: r12
REAL(8) :: r,couptd_local,width
INTEGER :: i,j,OMP_GET_NUM_THREADS
REAL(8), intent(out) :: couptd

! Description of the Variables:
! r12 is the distance between two points of the grid
! couptd is the coupling

couptd = 0.
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

couptd=couptd*dva*dvd

end Subroutine coultdc
