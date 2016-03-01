Subroutine diptrde(trden,dvx,dvy,dvz,origin,nx,ny,nz,mu)

    IMPLICIT None
! Input Variables
    INTEGER, intent(in) :: nx,ny,nz
	REAL(8), intent(in) :: dvx,dvy,dvz
	REAL(8), DIMENSION(3), intent(in) :: origin
	REAL(8), DIMENSION(nx,ny,nz), intent(in) :: trden
! Internal and output Variables
	INTEGER :: ix,iy,iz,nthrd,OMP_GET_NUM_THREADS
	REAL(8) :: dv,px,py,pz
	REAL(8), DIMENSION(nx,ny,nz) :: trd
    REAL(8) :: mux,muy,muz
    REAL(8), DIMENSION(3), intent(out) :: mu

dv=dvx*dvy*dvz
mux=0.d0
muy=0.d0
muz=0.d0

do iz=1,nz
	do iy=1,ny
		do ix=1,nx
	        px=origin(1)+(ix-1)*dvx
		    py=origin(2)+(iy-1)*dvy
			pz=origin(3)+(iz-1)*dvz
            trd(ix,iy,iz)=trden(ix,iy,iz)*dv
            mux=mux+px*trd(ix,iy,iz)
            muy=muy+py*trd(ix,iy,iz)
            muz=muz+pz*trd(ix,iy,iz)
        end do
    end do
end do

mu(1)=mux
mu(2)=muy
mu(3)=muz

end Subroutine diptrde


Subroutine couptr(trdena,dvxa,dvya,dvza,origina,trdend,dvxd,dvyd,dvzd,origind,nxa,nya,nza,nxd,nyd,nzd,thresh,couptd)

    IMPLICIT None
! Input Variables
    INTEGER, intent(in) :: nxa,nya,nza,nxd,nyd,nzd
	REAL(8), intent(in) :: dvxa,dvya,dvza,thresh
	REAL(8), intent(in) :: dvxd,dvyd,dvzd
	REAL(8), DIMENSION(3), intent(in) :: origina
	REAL(8), DIMENSION(3), intent(in) :: origind
	REAL(8), DIMENSION(nxa,nya,nza), intent(in) :: trdena
	REAL(8), DIMENSION(nxd,nyd,nzd), intent(in) :: trdend
! Internal and output Variables
    REAL(8) :: autown,start,finish
	INTEGER :: ixa,iya,iza,ixd,iyd,izd,nthrd,OMP_GET_NUM_THREADS
	REAL(8) :: dump,r12,dva,dvd,pxd,pyd,pzd,pxa,pya,pza
	REAL(8), DIMENSION(nxa,nya,nza) :: trda
	REAL(8), DIMENSION(nxd,nyd,nzd) :: trdd
    REAL(8), intent(out) :: couptd

! Description of the Variables:
! ixa,iya,iza,ixd,iyd,izd are the cycles indices 
! nx, ny, nz are the points of the grid
! dvx*dvy*dvz is the cubic volume
! origin is the starting point of the cube file
! dump is the dumping factor, it acts if the distance between two points is too short
! r12 is the distance between two points of the grid
! couptd is the coupling

couptd = 0.
autown = 2.194746d5

dva=dvxa*dvya*dvza
dvd=dvxd*dvyd*dvzd

! cycles along the six coordinates of the two density matrices

!$omp parallel do &
!$omp default(shared) &
!$omp private(ixa,iya,iza,pxa,pya,pza,ixd,iyd,izd,pxd,pyd,pzd,r12,dump)
do iza=1,nza
	pza=origina(3)+(iza-1)*dvza
	do iya=1,nya
		pya=origina(2)+(iya-1)*dvya
		do ixa=1,nxa
			if (abs(trdena(ixa,iya,iza)) < thresh) cycle
			pxa=origina(1)+(ixa-1)*dvxa
            !trda(ixa,iya,iza)=trdena(ixa,iya,iza)*dva
			do izd=1,nzd
				pzd=origind(3)+(izd-1)*dvzd
				do iyd=1,nyd
                    pyd=origind(2)+(iyd-1)*dvyd
					do ixd=1,nxd
						if (abs(trdend(ixd,iyd,izd)) < thresh) cycle
                        pxd=origind(1)+(ixd-1)*dvxd
						r12=(pxa-pxd)**2+(pya-pyd)**2+(pza-pzd)**2
						r12=dsqrt(r12)
                        dump=1.d0
						!if (r12.lt.5.d0) then
                        !    dump=exp(-(r12-5.d0)**2/1.d0)
                        !end if
                        !$omp critical
						couptd=couptd+(dump*trdena(ixa,iya,iza)*trdend(ixd,iyd,izd))/r12
                        !$omp end critical 
					end do
				end do
			end do
		end do
	end do
end do

!$omp end parallel do
couptd=couptd*dva*dvd*autown

end Subroutine couptr
