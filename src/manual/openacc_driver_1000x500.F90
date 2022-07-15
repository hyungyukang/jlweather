module openacc_driver
USE, INTRINSIC :: ISO_C_BINDING
USE OPENACC

INTEGER (C_INT64_T ), PARAMETER :: NX = 1000
INTEGER (C_INT64_T ), PARAMETER :: NZ = 500
REAL (C_DOUBLE ), PARAMETER :: DX = 20.0
REAL (C_DOUBLE ), PARAMETER :: DZ = 20.0
INTEGER (C_INT64_T ), PARAMETER :: HS = 2
INTEGER (C_INT64_T ), PARAMETER :: NUM_VARS = 4
REAL (C_DOUBLE ), PARAMETER :: C0 = 27.562941092972594
REAL (C_DOUBLE ), PARAMETER :: GAMMA = 1.400278940027894
REAL (C_DOUBLE ), PARAMETER :: P0 = 100000.0
REAL (C_DOUBLE ), PARAMETER :: HV_BETA = 0.05
REAL (C_DOUBLE ), PARAMETER :: GRAV = 9.8
REAL (C_DOUBLE ), PARAMETER :: RD = 287.0
REAL (C_DOUBLE ), PARAMETER :: CP = 1004.0
REAL (C_DOUBLE ), PARAMETER :: CV = 717.0
INTEGER (C_INT64_T ), PARAMETER :: ID_DENS = 1
INTEGER (C_INT64_T ), PARAMETER :: ID_UMOM = 2
INTEGER (C_INT64_T ), PARAMETER :: ID_WMOM = 3
INTEGER (C_INT64_T ), PARAMETER :: ID_RHOT = 4
INTEGER (C_INT64_T ), PARAMETER :: STEN_SIZE = 4

public jai_allocate, jai_updateto_init, jai_tend_x, jai_tend_y, jai_reductions, jai_deallocate 

contains

INTEGER (C_INT64_T) FUNCTION jai_allocate(state,statetmp,flux,tend,hy_dens_cell,hy_dens_theta_cell,hy_dens_int,hy_dens_theta_int,hy_pressure_int) BIND(C, name="jai_allocate")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state
REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: statetmp
REAL (C_DOUBLE ), DIMENSION(1:1001, 1:501, 1:4), INTENT(IN) :: flux
REAL (C_DOUBLE ), DIMENSION(1:1000, 1:500, 1:4), INTENT(IN) :: tend
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_cell
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_theta_cell
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_theta_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_pressure_int

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

!$acc enter data create(state)

!$acc enter data create(statetmp)

!$acc enter data create(flux)

!$acc enter data create(tend)

!$acc enter data create(hy_dens_cell)

!$acc enter data create(hy_dens_theta_cell)

!$acc enter data create(hy_dens_int)

!$acc enter data create(hy_dens_theta_int)

!$acc enter data create(hy_pressure_int)


jai_allocate = JAI_ERRORCODE

END FUNCTION

INTEGER (C_INT64_T) FUNCTION jai_updateto(state) BIND(C, name="jai_updateto")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

!$acc update device(state)

jai_updateto = JAI_ERRORCODE

END FUNCTION


INTEGER (C_INT64_T) FUNCTION jai_updateto_init(hy_dens_cell,hy_dens_theta_cell,hy_dens_int,hy_dens_theta_int,hy_pressure_int) BIND(C, name="jai_updateto_init")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_cell
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_theta_cell
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_theta_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_pressure_int

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0


!$acc update device(hy_dens_cell)

!$acc update device(hy_dens_theta_cell)

!$acc update device(hy_dens_int)

!$acc update device(hy_dens_theta_int)

!$acc update device(hy_pressure_int)

jai_updateto_init = JAI_ERRORCODE

END FUNCTION

INTEGER (C_INT64_T) FUNCTION jai_reductions(state,hy_dens_cell,hy_dens_theta_cell,glob) BIND(C, name="jai_reductions")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_cell
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_theta_cell
REAL (C_DOUBLE ), DIMENSION(1:2), INTENT(OUT) :: glob

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

    integer :: i, k, ierr
    real(C_DOUBLE) :: r,u,w,th,p,t,ke,ie, mass, te
    mass = 0.
    te = 0.
    !$acc parallel loop collapse(2) reduction(+:mass,te) present(state, hy_dens_cell, hy_dens_theta_cell)
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
    glob(1) = mass
    glob(2) = te

jai_reductions = JAI_ERRORCODE


END FUNCTION


INTEGER (C_INT64_T) FUNCTION jai_tend_x(state,dt,hy_dens_cell,hy_dens_theta_cell,flux,tend) BIND(C, name="jai_tend_x")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state
REAL (C_DOUBLE ), INTENT(IN) :: dt
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_cell
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_theta_cell
REAL (C_DOUBLE ), DIMENSION(1:1001, 1:501, 1:4), INTENT(OUT) :: flux
REAL (C_DOUBLE ), DIMENSION(1:1000, 1:500, 1:4), INTENT(OUT) :: tend

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

    integer :: i,k,ll,s
    real(C_DOUBLE) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dx / (16*dt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: THREAD ME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Compute fluxes in the x-direction for each cell
    !$acc parallel loop collapse(2) private(stencil,vals,d3_vals) &
    !$acc& present(state, hy_dens_cell, hy_dens_theta_cell, flux)
    do k = 1 , nz
      do i = 1 , nx+1
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i-hs-1+s,k,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state (for artificial viscosity)
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo
        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_cell(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_cell(k) ) / r
        p = C0*(r*t)**gamma
        !Compute the flux vector
        flux(i,k,ID_DENS) = r*u     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*u*u+p - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*u*w   - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*u*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: THREAD ME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Use the fluxes to compute tendencies for each cell
    !$acc parallel loop collapse(3) present(flux, tend)
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i+1,k,ll) - flux(i,k,ll) ) / dx
        enddo
      enddo
    enddo

jai_tend_x = JAI_ERRORCODE

END FUNCTION

INTEGER (C_INT64_T) FUNCTION jai_tend_z(state,dt,hy_dens_int,hy_dens_theta_int,hy_pressure_int,flux,tend) BIND(C, name="jai_tend_z")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state
REAL (C_DOUBLE ), INTENT(IN) :: dt
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_theta_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_pressure_int
REAL (C_DOUBLE ), DIMENSION(1:1001, 1:501, 1:4), INTENT(OUT) :: flux
REAL (C_DOUBLE ), DIMENSION(1:1000, 1:500, 1:4), INTENT(OUT) :: tend

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

    integer :: i,k,ll,s
    real(C_DOUBLE) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dz / (16*dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: THREAD ME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Compute fluxes in the x-direction for each cell
    !$acc parallel loop collapse(2) private(stencil,vals,d3_vals) &
    !$acc& present(state, hy_dens_int, hy_dens_theta_int, hy_pressure_int, flux)
    do k = 1 , nz+1
      do i = 1 , nx
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i,k-hs-1+s,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo
        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals(ID_DENS) + hy_dens_int(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_int(k) ) / r
        p = C0*(r*t)**gamma - hy_pressure_int(k)
        !Enforce vertical boundary condition and exact mass conservation
        if (k == 1 .or. k == nz+1) then
          w                = 0
          d3_vals(ID_DENS) = 0
        endif
        !Compute the flux vector with hyperviscosity
        flux(i,k,ID_DENS) = r*w     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*w*u   - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*w*w+p - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*w*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: THREAD ME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Use the fluxes to compute tendencies for each cell
    !$acc parallel loop collapse(3) present(flux, state, tend)
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i,k+1,ll) - flux(i,k,ll) ) / dz
          if (ll == ID_WMOM) then
            tend(i,k,ID_WMOM) = tend(i,k,ID_WMOM) - state(i,k,ID_DENS)*grav
          endif
        enddo
      enddo
    enddo

jai_tend_z = JAI_ERRORCODE


END FUNCTION

INTEGER (C_INT64_T) FUNCTION jai_updatefrom(tend) BIND(C, name="jai_updatefrom")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(1:1000, 1:500, 1:4), INTENT(IN) :: tend
INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

!$acc update host(tend)

jai_updatefrom = JAI_ERRORCODE

END FUNCTION

INTEGER (C_INT64_T) FUNCTION jai_deallocate(state,statetmp,flux,tend,hy_dens_cell,hy_dens_theta_cell,hy_dens_int,hy_dens_theta_int,hy_pressure_int) BIND(C, name="jai_deallocate")
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: state
REAL (C_DOUBLE ), DIMENSION(-1:1002, -1:502, 1:4), INTENT(IN) :: statetmp
REAL (C_DOUBLE ), DIMENSION(1:1001, 1:501, 1:4), INTENT(IN) :: flux
REAL (C_DOUBLE ), DIMENSION(1:1000, 1:500, 1:4), INTENT(IN) :: tend
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_cell
REAL (C_DOUBLE ), DIMENSION(-1:502), INTENT(IN) :: hy_dens_theta_cell
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_dens_theta_int
REAL (C_DOUBLE ), DIMENSION(1:501), INTENT(IN) :: hy_pressure_int

INTEGER (C_INT64_T) :: JAI_ERRORCODE  = 0

!$acc exit data delete(state)

!$acc exit data delete(statetmp)

!$acc exit data delete(flux)

!$acc exit data delete(tend)

!$acc exit data delete(hy_dens_cell)

!$acc exit data delete(hy_dens_theta_cell)

!$acc exit data delete(hy_dens_int)

!$acc exit data delete(hy_dens_theta_int)

!$acc exit data delete(hy_pressure_int)


jai_deallocate = JAI_ERRORCODE

END FUNCTION

end module


