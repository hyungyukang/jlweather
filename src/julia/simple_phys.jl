#-----------------------------------------------------------------------
#
#  Date:  April 26, 2016 (Version 6)
#
#  Simple Physics Package
#
#  SIMPLE_PHYSICS includes large-scale precipitation, surface fluxes and
#  boundary-leyer mixing. The processes are time-split in that order.
#  A partially implicit formulation is used to foster numerical
#  stability. The routine assumes that the model levels are ordered
#  in a top-down approach, e.g. level 1 denotes the uppermost full model
#  level.
#
#  This routine is based on an implementation which was developed for
#  the NCAR Community Atmosphere Model (CAM). Adjustments for other
#  models may be necessary.
#
#  The routine provides both updates of the state variables u, v, T, q
#  (these are local copies of u,v,T,q within this physics routine) and
#  also collects their time tendencies. The latter might be used to
#  couple the physics and dynamics in a process-split way. For a
#  time-split coupling, the final state should be given to the
#  dynamical core for the next time step.
#
# Test:      0 = Reed and Jablonowski (2011) tropical cyclone test
#            1 = Moist baroclinic instability test
# RJ2012_precip:
#         true  = Turn on Reed and Jablonowski (2012) precip scheme 
#         false = Turn off Reed and Jablonowski (2012) precip scheme
# TC_PBL_mod:
#         true  = Turn on George Bryan PBL mod for tropical cyclone test
#         false = Turn off George Bryan PBL mod [i.e., run as in Reed and Jablonowski (2012))
#
#  SUBROUTINE SIMPLE_PHYSICS(pcols, pver, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test, RJ2012_precip, TC_PBL_mod)
#
#  Input variables:
#     pcols  - number of atmospheric columns (#)
#     pver   - number of model levels (#)
#     dtime  - time step (s)
#     lat    - latitude (radians)
#     t      - temperature at model levels (K)
#     q      - specific humidity at model levels (gm/gm)
#     u      - zonal wind at model levels (m/s)
#     v      - meridional wind at model levels (m/s)
#     pmid   - pressure at model levels (Pa)
#     pint   - pressure at interfaces (Pa)
#     pdel   - layer thickness (Pa)
#     rpdel  - reciprocal of layer thickness (1/Pa)
#     ps     - surface pressure (Pa)
#     test   - test case to use for sea-surface temperatures
#     RJ2012_precip - RJ2012 precip flag
#     TC_PBL_mod    - PCL modification for TC test 
#
#  Output variables:
#     Increments are added into t, q, u, v, pmid, pint, pdel, rpdel and ps
#     which are returned to the routine from which SIMPLE_PHYSICS was
#     called.  Precpitation is returned via precl.
#
#  Change log:
#  v2: removal of some NCAR CAM-specific 'use' associations
#  v3: corrected precl[i] computation, the precipitation rate is now
#      computed via a vertical integral, the previous single-level
#      computation in v2 was a bug
#  v3: corrected dtdt[i,1) computation, the term '-[i,1)' was missing
#      the temperature variable: '-t[i,1)'
#  v4: modified and enhanced parameter list to make the routine truly
#      standalone, the number of columns and vertical levels have been
#      added: pcols, pver
#  v4: 'ncol' has been removed, 'pcols' is used instead
#  v5: the sea surface temperature (SST) field Tsurf is now an array,
#      the SST now depends on the latitude
#  v5: addition of the latitude array 'lat' and the flag 'test' in the
#      parameter list
#      if test = 0: constant SST is used, correct setting for the
#                   tropical cyclone test
#      if test = 1: newly added latitude-dependent SST is used,
#                   correct setting for the moist baroclinic wave test
#                   with simple-physics
#  v6: addition of flags for a modified PBL for the TC test and
#      to turn off large-scale condensation scheme when using Kessler physics
#      Included virtual temperature in density calculation in PBL scheme
#      Also, included the virtual temperature, instead of temperature, for
#      the calculation of rho in the PBL scheme
#      (v6_1) Minor specification and generalization fixes.
#      
# Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone
#            simulations of intermediate complexity: A test case for AGCMs,
#            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
#-----------------------------------------------------------------------

function simple_physics!(pcols::Int, 
                         pver::Int,
                         dtime::Float64,
                         t::Array{Float64,2},
                         q::Array{Float64,2},
                         u::Array{Float64,2},
                         v::Array{Float64,2},
                         pmid::Array{Float64,2},
                         pint::Array{Float64,2},
                         pdel::Array{Float64,2},
                         rpdel::Array{Float64},
                         ps::Vector{Float64},
                         precl::Vector{Float64})
#                        test,
#                        RJ2012_precip,
#                        TC_PBL_mod)

# ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
# ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

#  implicit none

#  integer, parameter :: r8 = selected_real_kind(12)

#
# Input arguments - MODEL DEPENDENT
#
#  integer, intent[in)  :: pcols        ! Set number of atmospheric columns
#  integer, intent[in)  :: pver         ! Set number of model levels
#  real(r8), intent[in) :: dtime        ! Set model physics timestep
#  integer, intent[in)  :: test         ! Test number
#  logical, intent[in)  :: RJ2012_precip
#  logical, intent[in)  :: TC_PBL_mod

#
# Input/Output arguments
#
#  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
#
#  real(r8), intent[inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
#  real(r8), intent[inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
#  real(r8), intent[inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
#  real(r8), intent[inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
#  real(r8), intent[inout) :: pmid(pcols,pver)   ! Pressure is full-model level (Pa)
#  real(r8), intent[inout) :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa)
#  real(r8), intent[inout) :: pdel(pcols,pver)   ! Layer thickness (Pa)
#  real(r8), intent[inout) :: rpdel(pcols,pver)  ! Reciprocal of layer thickness (1/Pa)
#  real(r8), intent[inout) :: ps(pcols)          ! Surface Pressue (Pa)

# Output arguments
#
#  real(r8), intent(out) :: precl(pcols)         ! Precipitation rate (m_water / s)

#
#---------------------------Local workspace-----------------------------
#

# Integers for loops

#  integer  i,k                         ! Longitude, level indices

# Physical Constants - Many of these may be model dependent

#  real(r8) gravit                      ! Gravity
#  real(r8) rair                        ! Gas constant for dry air
#  real(r8) cpair                       ! Specific heat of dry air
#  real(r8) latvap                      ! Latent heat of vaporization
#  real(r8) rh2o                        ! Gas constant for water vapor
#  real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
#  real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
#  real(r8) a                           ! Reference Earth's Radius (m)
#  real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
#  real(r8) pi                          ! pi

# Simple Physics Specific Constants 

#++++++++
#  real(r8) Tsurf(pcols)                ! Sea Surface Temperature (constant for tropical cyclone)
#++++++++                                 Tsurf needs to be dependent on latitude for the
#                                       ! moist baroclinic wave test, adjust

#  real(r8) SST_TC                      ! Sea Surface Temperature for tropical cyclone test
#  real(r8) T0                          ! Control temp for calculation of qsat
#  real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
#  real(r8) rhow                        ! Density of Liquid Water
#  real(r8) p0                          ! Constant for calculation of potential temperature
#  real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
#  real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
#  real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
#  real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
#  real(r8) C                           ! Drag coefficient for sensible heat and evaporation
#  real(r8) T00                         ! Horizontal mean T at surface for moist baro test
#  real(r8) u0                          ! Zonal wind constant for moist baro test
#  real(r8) latw                        ! halfwidth for  for baro test
#  real(r8) eta0                        ! Center of jets (hybrid) for baro test
#  real(r8) etav                        ! Auxiliary variable for baro test
#  real(r8) q0                          ! Maximum specific humidity for baro test
#  real(r8) kappa                       ! von Karman constant

# Physics Tendency Arrays
# real(r8) dtdt(pcols,pver)             ! Temperature tendency
# real(r8) dqdt(pcols,pver)             ! Specific humidity tendency
# real(r8) dudt(pcols,pver)             ! Zonal wind tendency
# real(r8) dvdt(pcols,pver)             ! Meridional wind tendency

  dtdt = Array{Float64}(undef, pcols, pver)
  dqdt = Array{Float64}(undef, pcols, pver)
  dudt = Array{Float64}(undef, pcols, pver)
  dvdt = Array{Float64}(undef, pcols, pver)

# Temporary variables for tendency calculations

#  real(r8) tmp                         ! Temporary
#  real(r8) qsat                        ! Saturation vapor pressure
#  real(r8) qsats                       ! Saturation vapor pressure of SST

# Variables for Boundary Layer Calculation

   wind = 

#  real(r8) wind(pcols)                 ! Magnitude of Wind
#  real(r8) Cd(pcols)                   ! Drag coefficient for momentum
#  real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
#  real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
#  real(r8) rho                         ! Density at lower/upper interface
#  real(r8) za(pcols)                   ! Heights at midpoints of first model level
#  real(r8) zi(pcols,pver+1)            ! Heights at model interfaces
#  real(r8) dlnpint                     ! Used for calculation of heights
#  real(r8) pbltop                      ! Top of boundary layer
#  real(r8) zpbltop                     ! Top of boundary layer for George Bryan Modifcation
#  real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
#  real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
#  real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
#  real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
#  real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
#  real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
#  real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
#  real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
#  real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
#  real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
#  real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme

   wind = Vector{Float64}(undef, pcols)       # Magnitude of Wind
   Cd = Vector{Float64}(undef, pcols)         # Drag coefficient for momentum
   Km = Array{Float64}(undef, pcols,pver+1)   # Eddy diffusivity for boundary layer calculations 
   Ke = Array{Float64}(undef, pcols,pver+1)   # Eddy diffusivity for boundary layer calculations
#  rho                                        # Density at lower/upper interface
   za = Vector{Float64}(undef, pcols)         # Heights at midpoints of first model level
   zi = Array{Float64}(undef, pcols,pver+1)   # Heights at model interfaces
#  dlnpint                                    # Used for calculation of heights
#  pbltop                                     # Top of boundary layer
#  zpbltop                                    # Top of boundary layer for George Bryan Modifcation
#  pblconst                                   # Constant for the calculation of the decay of diffusivity 
   CA = Array{Float64}(undef, pcols,pver)     # Matrix Coefficents for PBL Scheme 
   CC = Array{Float64}(undef, pcols,pver)     # Matrix Coefficents for PBL Scheme 
   CE = Array{Float64}(undef, pcols,pver+1)   # Matrix Coefficents for PBL Scheme
   CAm = Array{Float64}(undef, pcols,pver)    # Matrix Coefficents for PBL Scheme
   CCm = Array{Float64}(undef, pcols,pver)    # Matrix Coefficents for PBL Scheme
   CEm = Array{Float64}(undef, pcols,pver+1)  # Matrix Coefficents for PBL Scheme
   CFu = Array{Float64}(undef, pcols,pver+1)  # Matrix Coefficents for PBL Scheme
   CFv = Array{Float64}(undef, pcols,pver+1)  # Matrix Coefficents for PBL Scheme
   CFt = Array{Float64}(undef, pcols,pver+1)  # Matrix Coefficents for PBL Scheme
   CFq = Array{Float64}(undef, pcols,pver+1)  # Matrix Coefficents for PBL Scheme


# Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
# conserve the mass of the dry air

   qini = Array{Float64}(undef, pcols,pver)   # Initial specific humidity

#-----------------------------------------------------------------------------
#
#Physical Constants - MAY BE MODEL DEPENDENT
#
#-----------------------------------------------------------------------------

   gravit = 9.80616                  # Gravity (9.80616 m/s^2)
   rair   = 287.0                    # Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3                 # Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6                    # Latent heat of vaporization (J/kg)
   rh2o   = 461.5                    # Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                # Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1.0        # Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0                # Reference Earth's Radius (m)
   omega  = 7.29212e-5               # Reference rotation rate of the Earth (s^-1)
#  pi     = 4.0*atan(1.0)            # pi

#-----------------------------------------------------------------------------
#
# Local Constants for Simple Physics
#
#-----------------------------------------------------------------------------
      C        = 0.0011      ! From Simth and Vogl 2008
      SST_TC   = 302.15      ! Constant Value for SST
      T0       = 273.16      ! control temp for calculation of qsat
      e0       = 610.78      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0      ! Density of Liquid Water 
      Cd0      = 0.0007      ! Constant for Cd calc. Simth and Vogl 2008
      Cd1      = 0.000065    ! Constant for Cd calc. Simth and Vogl 2008
      Cm       = 0.002       ! Constant for Cd calc. Simth and Vogl 2008
      v20      = 20.0        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0    ! Constant for potential temp calculation
      pbltop   = 85000.0     ! Top of boundary layer in p
      zpbltop  = 1000.0      ! Top of boundary layer in z
      pblconst = 10000.0     ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0       ! Horizontal mean T at surface for moist baro test
      u0       = 35.0        ! Zonal wind constant for moist baro test
      latw     = 2.0*pi/9.0 ! Halfwidth for  for baro test
      eta0     = 0.252       ! Center of jets (hybrid) for baro test
      etav     = (1.0-eta0)*0.5*pi ! Auxiliary variable for baro test
      q0       = 0.021       ! Maximum specific humidity for baro test
      kappa    = 0.4         ! von Karman constant

#-----------------------------------------------------------------------------
#
# Definition of local arrays
#
#-----------------------------------------------------------------------------
#
# Calculate hydrostatic height
#
     for i in 1:pcols
        dlnpint = log(ps[i]) - log(pint[i,pver])  # ps[i] is identical to pint[i,pver+1), note: this is the correct sign (corrects typo in JAMES paper) 
        za[i] = rair/gravit*t[i,pver]*(1.0+zvir*q[i,pver])*0.50*dlnpint
        zi[i,pver+1] = 0.00
     end
#
# Set Initial Specific Humidity
#
     qini = q
#
# Set Sea Surface Temperature (constant for tropical cyclone)
# Tsurf needs to be dependent on latitude for moist baroclinic wave test
# Tsurf needs to be constant for tropical cyclone test
#
     if (test == 1) # Moist Baroclinic Wave Test
        for i in 1:pcols
           Tsurf[i] = (T00 + pi*u0/rair * 1.50 * sin(etav) * (cos(etav))**0.50 *                 &
                     ((-2.0*(sin(lat[i]))**6 * ((cos(lat[i]))**2 + 1.0/3.0) + 10.0/63.0)* &
                     u0 * (cos(etav))**1.50  +                                                    &
                     (8.0/5.0*(cos(lat[i]))**3 * ((sin(lat[i]))**2 + 2.0/3.0) - pi/4.0)*a*omega*0.50 ))/ &
                     (1.0+zvir*q0*exp(-(lat[i]/latw)**4))
        end
     elseif (test == 0) # Tropical Cyclone Test
        for i in 1:pcols
           Tsurf[i] = SST_TC
        end
     end

#-----------------------------------------------------------------------------
#
# Set initial physics time tendencies and precipitation field to zero
#
#-----------------------------------------------------------------------------
     dtdt = 0.0            # initialize temperature tendency with zero
     dqdt = 0.0            # initialize specific humidity tendency with zero
     dudt = 0.0            # initialize zonal wind tendency with zero
     dvdt = 0.0            # initialize meridional wind tendency with zero
     precl = 0.0           # initialize precipitation rate with zero

#-----------------------------------------------------------------------------
#
# Large-Scale Condensation and Precipitation
#
#-----------------------------------------------------------------------------

      if (RJ2012_precip)
#
# Calculate Tendencies
#
      for k in 1:pver
         for i in 1:pcols
            qsat = epsilo*e0/pmid[i,k]*exp(-latvap/rh2o*((1.0/t[i,k])-1.0/T0))
            if (q[i,k] > qsat)
               tmp  = 1.0/dtime*(q[i,k]-qsat)/(1.0+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t[i,k]**2)))
               dtdt[i,k] = dtdt[i,k]+latvap/cpair*tmp
               dqdt[i,k] = dqdt[i,k]-tmp
               precl[i]  = precl[i]+tmp*pdel[i,k]/(gravit*rhow)
            end
         end
      end
#
# Update moisture and temperature fields from Larger-Scale Precipitation Scheme
#
      for k in 1:pver
         for i in 1:pcols
            t[i,k] =  t[i,k] + dtdt[i,k]*dtime
            q[i,k] =  q[i,k] + dqdt[i,k]*dtime
         end
      end

#-----------------------------------------------------------------------------
# Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
#
# note: The variables, as done in many GCMs, are written to the history file
#       after the moist physics process.  This ensures that the moisture fields
#       are somewhat in equilibrium.
#-----------------------------------------------------------------------------
  #  call diag_phys_writeout(state)   ! This is for CESM/CAM
    
      end
     
#-----------------------------------------------------------------------------
#
# Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
# sensible heat and latent heat
#
# We are using Simplified Ekman theory to compute the diffusion coefficients
# Kx for the boundary-layer mixing. The Kx values are calculated at each time step
# and in each column.
#
#-----------------------------------------------------------------------------
# Compute magnitude of the wind and drag coeffcients for turbulence scheme
#
     for i in 1:pcols
        wind[i] = sqrt(u[i,pver]**2+v[i,pver]**2)
     end
     for i in 1:pcols
        if( wind[i] .lt. v20)
           Cd[i] = Cd0+Cd1*wind[i] 
        else
           Cd[i] = Cm
        end
     end

     if (TC_PBL_mod) #Bryan TC PBL Modification 
     for k in pver:-1:1
        for i in 1:pcols
           dlnpint = log(pint[i,k+1]) - log(pint[i,k])
           zi[i,k] = zi[i,k+1]+rair/gravit*t[i,k]*(1.0+zvir*q[i,k])*dlnpint
           if( zi[i,k] <= zpbltop)
              Km[i,k] = kappa*sqrt(Cd[i])*wind[i]*zi[i,k]*(1.0-zi[i,k]/zpbltop)*(1.0-zi[i,k]/zpbltop)
              Ke[i,k] = kappa*sqrt(C)*wind[i]*zi[i,k]*(1.0-zi[i,k]/zpbltop)*(1.0-zi[i,k]/zpbltop) 
           else
              Km[i,k] = 0.00
              Ke[i,k] = 0.00
           end
        end
     end
     else # Reed and Jablonowski (2012) Configuration
     for k in 1:pver
        for i in 1:pcols
           if( pint[i,k] >= pbltop)
              Km[i,k] = Cd[i]*wind[i]*za[i] 
              Ke[i,k] = C*wind[i]*za[i]
           else
              Km[i,k] = Cd[i]*wind[i]*za[i]*exp(-(pbltop-pint[i,k])**2/(pblconst)**2)
              Ke[i,k] = C*wind[i]*za[i]*exp(-(pbltop-pint[i,k])**2/(pblconst)**2)
           end
        end
     end
     end
#-----------------------------------------------------------------------------
# Update the state variables u, v, t, q with the surface fluxes at the
# lowest model level, this is done with an implicit approach
# see Reed and Jablonowski (JAMES, 2012)
#
# Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
# Tsurf needs to be dependent on latitude for the
# moist baroclinic wave test 
#-----------------------------------------------------------------------------
     for i in 1:pcols
        qsats = epsilo*e0/ps[i]*exp(-latvap/rh2o*((1.0/Tsurf[i])-1.0/T0))
        dudt[i,pver] = dudt[i,pver] + (u[i,pver]/(1.0+Cd[i]*wind[i]*dtime/za[i])-u[i,pver])/dtime
        dvdt[i,pver] = dvdt[i,pver] + (v[i,pver]/(1.0+Cd[i]*wind[i]*dtime/za[i])-v[i,pver])/dtime
        u[i,pver] = u[i,pver]/(1.0+Cd[i]*wind[i]*dtime/za[i])
        v[i,pver] = v[i,pver]/(1.0+Cd[i]*wind[i]*dtime/za[i])
        dtdt[i,pver] = dtdt[i,pver] +((t[i,pver]+C*wind[i]*Tsurf[i]*dtime/za[i])/(1.0+C*wind[i]*dtime/za[i])-t[i,pver])/dtime 
        t[i,pver] = (t[i,pver]+C*wind[i]*Tsurf[i]*dtime/za[i])/(1.0+C*wind[i]*dtime/za[i])  
        dqdt[i,pver] = dqdt[i,pver] +((q[i,pver]+C*wind[i]*qsats*dtime/za[i])/(1.0+C*wind[i]*dtime/za[i])-q[i,pver])/dtime
        q[i,pver] = (q[i,pver]+C*wind[i]*qsats*dtime/za[i])/(1.0+C*wind[i]*dtime/za[i])
     end

#-----------------------------------------------------------------------------
# Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
#-----------------------------------------------------------------------------
# Calculate Diagonal Variables for Implicit PBL Scheme
#
      for k in 1:pver-1
         for i in 1:pcols
            rho = (pint[i,k+1]/(rair*(t[i,k+1]*(1.0+zvir*q[i,k+1])+t[i,k]*(1.0+zvir*q[i,k]))/2.00)) 
            CAm[i,k] = rpdel[i,k]*dtime*gravit*gravit*Km[i,k+1]*rho*rho/(pmid[i,k+1]-pmid[i,k])    
            CCm[i,k+1] = rpdel[i,k+1]*dtime*gravit*gravit*Km[i,k+1]*rho*rho/(pmid[i,k+1]-pmid[i,k])
            CA[i,k] = rpdel[i,k]*dtime*gravit*gravit*Ke[i,k+1]*rho*rho/(pmid[i,k+1]-pmid[i,k])
            CC[i,k+1] = rpdel[i,k+1]*dtime*gravit*gravit*Ke[i,k+1]*rho*rho/(pmid[i,k+1]-pmid[i,k])
         end
      end
      for i in 1:pcols
         CAm[i,pver] = 0.0
         CCm[i,1] = 0.0
         CEm[i,pver+1] = 0.0
         CA[i,pver] = 0.0
         CC[i,1] = 0.0
         CE[i,pver+1] = 0.0
         CFu[i,pver+1] = 0.0
         CFv[i,pver+1] = 0.0
         CFt[i,pver+1] = 0.0
         CFq[i,pver+1] = 0.0 
      end
      for i in 1:pcols
         for k in pver:-1:1
            CE[i,k] = CC[i,k]/(1.0+CA[i,k]+CC[i,k]-CA[i,k]*CE[i,k+1]) 
            CEm[i,k] = CCm[i,k]/(1.0+CAm[i,k]+CCm[i,k]-CAm[i,k]*CEm[i,k+1])
            CFu[i,k] = (u[i,k]+CAm[i,k]*CFu[i,k+1])/(1.0+CAm[i,k]+CCm[i,k]-CAm[i,k]*CEm[i,k+1])
            CFv[i,k] = (v[i,k]+CAm[i,k]*CFv[i,k+1])/(1.0+CAm[i,k]+CCm[i,k]-CAm[i,k]*CEm[i,k+1])
            CFt[i,k] = ((p0/pmid[i,k])**(rair/cpair)*t[i,k]+CA[i,k]*CFt[i,k+1])/(1.0+CA[i,k]+CC[i,k]-CA[i,k]*CE[i,k+1]) 
            CFq[i,k] = (q[i,k]+CA[i,k]*CFq[i,k+1])/(1.0+CA[i,k]+CC[i,k]-CA[i,k]*CE[i,k+1])
       end
      end
#
# Calculate the updated temperaure and specific humidity and wind tendencies
#
# First we need to calculate the tendencies at the top model level
#
      for i in 1:pcols
            dudt[i,1]  = dudt[i,1]+(CFu[i,1]-u[i,1])/dtime
            dvdt[i,1]  = dvdt[i,1]+(CFv[i,1]-v[i,1])/dtime
            u[i,1]    = CFu[i,1]
            v[i,1]    = CFv[i,1]
            dtdt[i,1]  = dtdt[i,1]+(CFt[i,1]*(pmid[i,1]/p0)**(rair/cpair)-t[i,1])/dtime
            t[i,1]    = CFt[i,1]*(pmid[i,1]/p0)**(rair/cpair)
            dqdt[i,1]  = dqdt[i,1]+(CFq[i,1]-q[i,1])/dtime
            q[i,1]  = CFq[i,1]
      end

      for i in 1:pcols
         for k in 2:pver
            dudt[i,k]  = dudt[i,k]+(CEm[i,k]*u[i,k-1]+CFu[i,k]-u[i,k])/dtime
            dvdt[i,k]  = dvdt[i,k]+(CEm[i,k]*v[i,k-1]+CFv[i,k]-v[i,k])/dtime
            u[i,k]    = CEm[i,k]*u[i,k-1]+CFu[i,k] 
            v[i,k]    = CEm[i,k]*v[i,k-1]+CFv[i,k]
            dtdt[i,k]  = dtdt[i,k]+((CE[i,k]*t[i,k-1]*(p0/pmid[i,k-1])**(rair/cpair)+CFt[i,k])*(pmid[i,k]/p0)**(rair/cpair)-t[i,k])/dtime 
            t[i,k]    = (CE[i,k]*t[i,k-1]*(p0/pmid[i,k-1])**(rair/cpair)+CFt[i,k])*(pmid[i,k]/p0)**(rair/cpair)
            dqdt[i,k]  = dqdt[i,k]+(CE[i,k]*q[i,k-1]+CFq[i,k]-q[i,k])/dtime
            q[i,k]  = CE[i,k]*q[i,k-1]+CFq[i,k]
         end
      end

#-----------------------------------------------------------------------------
# Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC
#
# note: Care needs to be taken to ensure that the model conserves the total
#       dry air mass. Add your own routine here.
#-----------------------------------------------------------------------------
  #  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   return
end
