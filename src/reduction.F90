module reducemod
USE, INTRINSIC :: ISO_C_BINDING

public launch

contains

!INTEGER (C_INT64_T) FUNCTION launch(nx, nz, dx, dz, hs, NUM_VARS, &
!                                state) BIND(C, name="launch")
!
!    USE, INTRINSIC :: ISO_C_BINDING
!
!    implicit none
!
!    INTEGER (C_INT64_T), INTENT(IN) :: nx, nz, hs, NUM_VARS
!    REAL(C_DOUBLE), INTENT(IN) :: dx, dz
!    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS) :: state
!
!    INTEGER (C_INT64_T) :: res
!
!    print *, "Inside of reduction kernel"
!    print *, "scalars: ", nx, nz, dx, dz, hs, NUM_VARS
!
!    print *, "state(0, 0, 1) = ", state(0, 0, 1)
!
!    state(0, 0, 1) = 200.0
!
!    res = 0
!
!    launch = 0
!
!END FUNCTION

INTEGER (C_INT64_T) FUNCTION launch(nx, nz, dx, dz, hs, NUM_VARS, state, &
                   hy_dens_cell, hy_dens_theta_cell, C0, gamma, &
                   p0, rd, cp, cv, ID_DENS, ID_UMOM, ID_WMOM, &
                   ID_RHOT, glob) BIND(C, name="launch")

    USE, INTRINSIC :: ISO_C_BINDING

    implicit none

    INTEGER (C_INT64_T), INTENT(IN) :: nx, nz, hs, NUM_VARS, ID_DENS, ID_UMOM, ID_WMOM, ID_RHOT
    REAL(C_DOUBLE), INTENT(IN) :: dx, dz, C0, gamma,  p0, rd, cp, cv
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(1-hs:nz+hs) :: hy_dens_cell, hy_dens_theta_cell
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS) :: state
    REAL(C_DOUBLE), INTENT(OUT), DIMENSION(2) :: glob

    INTEGER (C_INT64_T) :: res, i, k, ierr

    REAL(C_DOUBLE) :: r,u,w,th,p,t,ke,ie
    REAL(C_DOUBLE) :: mass, te

    print *, "Inside of reduction kernel"

    res = 0

    mass = 0
    te   = 0

    do k = 1 , nz
      do i = 1 , nx
        r  =   state(i,k,ID_DENS) + hy_dens_cell(k)             ! Density
        u  =   state(i,k,ID_UMOM) / r                           ! U-wind
        w  =   state(i,k,ID_WMOM) / r                           ! W-wind
        th = ( state(i,k,ID_RHOT) + hy_dens_theta_cell(k) ) / r ! Potential Temperature (theta)
        p  = C0*(r*th)**gamma      ! Pressure
        t  = th / (p0/p)**(rd/cp)  ! Temperature
        ke = r*(u*u+w*w)           ! Kinetic Energy
        ie = r*cv*t                ! Internal Energy
        mass = mass + r            *dx*dz ! Accumulate domain mass
        te   = te   + (ke + r*cv*t)*dx*dz ! Accumulate domain total energy
      enddo
    enddo

    print *, "MASS K : ", mass, te

    glob(1) = mass
    glob(2) = te

    launch = res

END FUNCTION

end module
