using OffsetArrays
using Match
using MPI
using Debugger
using Printf

##############
# constants
##############
    
# julia command to link MPI.jl to system MPI installation
# julia -e 'ENV["JULIA_MPI_BINARY"]="system"; ENV["JULIA_MPI_PATH"]="/Users/8yk/opt/usr/local"; using Pkg; Pkg.build("MPI"; verbose=true)'
MPI.Init()
const COMM   = MPI.COMM_WORLD
const NRANKS = MPI.Comm_size(COMM)
const MYRANK = MPI.Comm_rank(COMM)

if length(ARGS) >= 4
    const SIM_TIME    = parse(Float64, ARGS[1])
    const NX_GLOB     = parse(Int64, ARGS[2])
    const NZ_GLOB     = parse(Int64, ARGS[3])
    const DATA_SPEC   = parse(Int64, ARGS[4])
else
    const SIM_TIME    = Float64(100.0)
    const NX_GLOB     = Int64(200)
    const NZ_GLOB     = Int64(100)
    const DATA_SPEC   = Int64(1)
end
    
const NPER  = Float64(NX_GLOB)/NRANKS
const I_BEG = trunc(Int, round(NPER* MYRANK)+1)
const I_END = trunc(Int, round(NPER*(MYRANK+1)))
const NX    = I_END - I_BEG + 1
const NZ = NZ_GLOB

const LEFT_RANK = MYRANK-1 == -1 ? NRANKS - 1 : MYRANK - 1 
const RIGHT_RANK = MYRANK+1 == NRANKS ? 0 : MYRANK + 1 

#Vertical direction isn't MPI-ized, so the rank's local values = the global values
K_BEG      = 1
MASTERPROC = (MYRANK == 0)

const HS          = 4
const NUM_VARS    = 4
const XLEN        = Float64(2.E4) # Length of the domain in the x-direction (meters)
const ZLEN        = Float64(1.E4) # Length of the domain in the z-direction (meters)
const CFL         = Float64(1.5)  # "Courant, Friedrichs, Lewy" number (for numerical stability)
const MAX_SPEED   = Float64(450.0)# Assumed maximum wave speed during the simulation (speed of sound + speed of wind) (meter / sec)
const DX          = XLEN / NX_GLOB
const DZ          = ZLEN / NZ_GLOB
const DT          = min(DX,DZ) / MAX_SPEED * CFL
const NQPOINTS    = 3
const CP          = Float64(1004.0) # Specific heat of dry air at constant pressure
const CV          = Float64(717.0)  # Specific heat of dry air at constant volume
const RD          = Float64(287.0)  # Dry air constant for equation of state (P=rho*rd*T)
const P0          = Float64(1.0E5)  # Standard pressure at the surface in Pascals
const C0          = Float64(27.5629410929725921310572974482)
const GAMMA       = Float64(1.40027894002789400278940027894)

const ID_DENS     = 1
const ID_UMOM     = 2
const ID_WMOM     = 3
const ID_RHOT     = 4
                    
const DIR_X       = 1 #Integer constant to express that this operation is in the x-direction
const DIR_Z       = 2 #Integer constant to express that this operation is in the z-direction

const DATA_SPEC_COLLISION       = 1
const DATA_SPEC_THERMAL         = 2
const DATA_SPEC_MOUNTAIN        = 3
const DATA_SPEC_TURBULENCE      = 4
const DATA_SPEC_DENSITY_CURRENT = 5
const DATA_SPEC_INJECTION       = 6

const qpoints     = Array{Float64}([0.112701665379258311482073460022E0 , 0.500000000000000000000000000000E0 , 0.887298334620741688517926539980E0])
const qweights    = Array{Float64}([0.277777777777777777777777777779E0 , 0.444444444444444444444444444444E0 , 0.277777777777777777777777777779E0])

##############
# functions
##############

"""
    main()

simulate weather-like flows.

# Examples
```julia-repl
julia> main()
```
"""
function main(args::Vector{String})

    ######################
    # top-level variables
    ######################

    #@show(args)
    
    etime = Float64(0.0)
    dt = DT
    
    state, statetmp, flux, tend, hy_dens_cell, hy_dens_theta_cell = init!()

    mass0, te0 = reductions(state, hy_dens_cell, hy_dens_theta_cell)

    output()

    # main loop
    elapsedtime = @elapsed while etime < SIM_TIME

        if etime + dt > SIM_TIME
            dt = SIM_TIME - etime
        end

        timestep!(state, statetmp, flux, tend, dt)

        etime = etime + dt

    end
    
    mass, te = reductions(state, hy_dens_cell, hy_dens_theta_cell)

    if MASTERPROC
        println( "CPU Time: $elapsedtime")
        @printf("d_mass: %f\n", (mass - mass0)/mass0)
        @printf("d_te:   %f\n", (te - te0)/te0)
    end
        
    finalize!(state)

end

function init!()
    
    if MASTERPROC
        println("nx_glob, nz_glob: $NX_GLOB $NZ_GLOB")
        println("dx, dz: $DX $DZ")
    end
        
    #println("nx, nz at $MYRANK: $NX($I_BEG:$I_END) $NZ($K_BEG:$NZ)")
    
    MPI.Barrier(COMM)
    
    _state      = zeros(Float64, NX+2*HS, NZ+2*HS, NUM_VARS) 
    state       = OffsetArray(_state, 1-HS:NX+HS, 1-HS:NZ+HS, 1:NUM_VARS)
    _statetmp   = Array{Float64}(undef, NX+2*HS, NZ+2*HS, NUM_VARS) 
    statetmp    = OffsetArray(_state, 1-HS:NX+HS, 1-HS:NZ+HS, 1:NUM_VARS)
    
    flux        = Array{Float64}(undef, NX+1, NZ+1, NUM_VARS) 
    tend        = Array{Float64}(undef, NX, NZ, NUM_VARS) 

    _hy_dens_cell = zeros(Float64, NZ+2*HS) 
    hy_dens_cell  = OffsetArray(_hy_dens_cell, 1-HS:NZ+HS)
    _hy_dens_theta_cell = zeros(Float64, NZ+2*HS) 
    hy_dens_theta_cell  = OffsetArray(_hy_dens_theta_cell, 1-HS:NZ+HS)   
    
    hy_dens_int         = Array{Float64}(undef, NZ+1)
    hy_dens_theta_int   = Array{Float64}(undef, NZ+1)
    hy_pressure_int     = Array{Float64}(undef, NZ+1)
    
    sendbuf_l   = Array{Float64}(undef, HS, NZ, NUM_VARS)
    sendbuf_r   = Array{Float64}(undef, HS, NZ, NUM_VARS)
    recvbuf_l   = Array{Float64}(undef, HS, NZ, NUM_VARS)
    recvbuf_r   = Array{Float64}(undef, HS, NZ, NUM_VARS)
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Initialize the cell-averaged fluid state via Gauss-Legendre quadrature
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for k in 1-HS:NZ+HS
      for i in 1-HS:NX+HS
        #Initialize the state to zero
        #Use Gauss-Legendre quadrature to initialize a hydrostatic balance + temperature perturbation
        for kk in 1:NQPOINTS
          for ii in 1:NQPOINTS
            #Compute the x,z location within the global domain based on cell and quadrature index
            x = (I_BEG-1 + i-0.5) * DX + (qpoints[ii]-0.5)*DX
            z = (K_BEG-1 + k-0.5) * DZ + (qpoints[kk]-0.5)*DZ

            #Set the fluid state based on the user's specification
            r, u, w, t, hr, ht = @match data_spec_int begin
                DATA_SPEC_COLLISION       => collision!(x,z)
                DATA_SPEC_THERMAL         => thermal!(x,z)
                DATA_SPEC_MOUNTAIN        => mountain_waves!(x,z)
                DATA_SPEC_TURBULENCE      => turbulence!(x,z)
                DATA_SPEC_DENSITY_CURRENT => density_current!(x,z)
                DATA_SPEC_INJECTION       => injection!(x,z)
            end

            #Store into the fluid state array
            state[i,k,ID_DENS] = state[i,k,ID_DENS] + r                         * qweights[ii]*qweights[kk]
            state[i,k,ID_UMOM] = state[i,k,ID_UMOM] + (r+hr)*u                  * qweights[ii]*qweights[kk]
            state[i,k,ID_WMOM] = state[i,k,ID_WMOM] + (r+hr)*w                  * qweights[ii]*qweights[kk]
            state[i,k,ID_RHOT] = state[i,k,ID_RHOT] + ( (r+hr)*(t+ht) - hr*ht ) * qweights[ii]*qweights[kk]
          end
        end
        for ll in 1:NUM_VARS
          statetmp[i,k,ll] = state[i,k,ll]
        end
      end
    end

    for k in 1-HS:NZ+HS
        for kk in 1:NQPOINTS
            z = (K_BEG-1 + k-0.5) * DZ + (qpoints[kk]-0.5)*DZ
            
            #Set the fluid state based on the user's specification
            r, u, w, t, hr, ht = @match data_spec_int begin
                DATA_SPEC_COLLISION       => collision!(0.0,z)
                DATA_SPEC_THERMAL         => thermal!(0.0,z)
                DATA_SPEC_MOUNTAIN        => mountain_waves!(0.0,z)
                DATA_SPEC_TURBULENCE      => turbulence!(0.0,z)
                DATA_SPEC_DENSITY_CURRENT => density_current!(0.0,z)
                DATA_SPEC_INJECTION       => injection!(0.0,z)
            end           

            hy_dens_cell[k]       = hy_dens_cell[k]       + hr    * qweights[kk]
            hy_dens_theta_cell[k] = hy_dens_theta_cell[k] + hr*ht * qweights[kk]
        end
    end
    
    #Compute the hydrostatic background state at vertical cell interfaces
    for k in 1:NZ+1
        z = (K_BEG-1 + k-1) * DZ
        #Set the fluid state based on the user's specification
        r, u, w, t, hr, ht = @match data_spec_int begin
            DATA_SPEC_COLLISION       => collision!(0.0,z)
            DATA_SPEC_THERMAL         => thermal!(0.0,z)
            DATA_SPEC_MOUNTAIN        => mountain_waves!(0.0,z)
            DATA_SPEC_TURBULENCE      => turbulence!(0.0,z)
            DATA_SPEC_DENSITY_CURRENT => density_current!(0.0,z)
            DATA_SPEC_INJECTION       => injection!(0.0,z)
        end                  

      hy_dens_int[k] = hr
      hy_dens_theta_int[k] = hr*ht
      hy_pressure_int[k] = C0*(hr*ht)^GAMMA
    end
                                                                
    @bp
    return state, statetmp, flux, tend, hy_dens_cell, hy_dens_theta_cell
end

function injection!(x::Float64, z::Float64)
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

function density_current!(x::Float64, z::Float64)
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

function turbulence!(x::Float64, z::Float64)
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

function mountain_waves!(x::Float64, z::Float64)
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

function thermal!(x::Float64, z::Float64)
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

function collision!(x::Float64, z::Float64)
    
    r  = Float64(0.0)
    t  = Float64(1.0)
    u  = Float64(0.0)
    w  = Float64(0.0)
    hr = Float64(2.0)
    ht = Float64(3.0)
    
    return r, u, w, t, hr, ht
end

#Performs a single dimensionally split time step using a simple low-storate three-stage Runge-Kutta time integrator
#The dimensional splitting is a second-order-accurate alternating Strang splitting in which the
#order of directions is alternated each time step.
#The Runge-Kutta method used here is defined as follows:
# q*     = q[n] + dt/3 * rhs(q[n])
# q**    = q[n] + dt/2 * rhs(q*  )
# q[n+1] = q[n] + dt/1 * rhs(q** )
function timestep!(state::OffsetArray{Float64, 3, Array{Float64, 3}},
                   statetmp::OffsetArray{Float64, 3, Array{Float64, 3}},
                   flux::Array{Float64, 3},
                   tend::Array{Float64, 3},
                   dt::Float64)
    
    local direction_switch = true
    
    if direction_switch
        
        #x-direction first
        semi_discrete_step!( state , state    , statetmp , dt / 3 , DIR_X , flux , tend )
        semi_discrete_step!( state , statetmp , statetmp , dt / 2 , DIR_X , flux , tend )
        semi_discrete_step!( state , statetmp , state    , dt / 1 , DIR_X , flux , tend )
        
        #z-direction second
        semi_discrete_step!( state , state    , statetmp , dt / 3 , DIR_Z , flux , tend )
        semi_discrete_step!( state , statetmp , statetmp , dt / 2 , DIR_Z , flux , tend )
        semi_discrete_step!( state , statetmp , state    , dt / 1 , DIR_Z , flux , tend )
    else
        
        #z-direction second
        semi_discrete_step!( state , state    , statetmp , dt / 3 , DIR_Z , flux , tend )
        semi_discrete_step!( state , statetmp , statetmp , dt / 2 , DIR_Z , flux , tend )
        semi_discrete_step!( state , statetmp , state    , dt / 1 , DIR_Z , flux , tend )
        
        #x-direction first
        semi_discrete_step!( state , state    , statetmp , dt / 3 , DIR_X , flux , tend )
        semi_discrete_step!( state , statetmp , statetmp , dt / 2 , DIR_X , flux , tend )
        semi_discrete_step!( state , statetmp , state    , dt / 1 , DIR_X , flux , tend )
    end
end

        
#Perform a single semi-discretized step in time with the form:
#state_out = state_init + dt * rhs(state_forcing)
#Meaning the step starts from state_init, computes the rhs using state_forcing, and stores the result in state_out
function semi_discrete_step!(stateinit::OffsetArray{Float64, 3, Array{Float64, 3}},
                   stateforcing::OffsetArray{Float64, 3, Array{Float64, 3}},
                   stateout::OffsetArray{Float64, 3, Array{Float64, 3}},
                   dt::Float64,
                   dir::Int,
                   flux::Array{Float64, 3},
                   tend::Array{Float64, 3})

    if dir == DIR_X
      #Set the halo values for this MPI task's fluid state in the x-direction
      set_halo_values_x!(stateforcing)
        
      #Compute the time tendencies for the fluid state in the x-direction
      compute_tendencies_x!(stateforcing,flux,tend)
    elseif dir == DIR_Z
      #Set the halo values for this MPI task's fluid state in the z-direction
      set_halo_values_z!(stateforcing)
        
      #Compute the time tendencies for the fluid state in the z-direction
      compute_tendencies_z!(stateforcing,flux,tend)
    end

    #Apply the tendencies to the fluid state
    for ll in 1:NUM_VARS
        for k in 1:NZ
            for i in 1:NX
                stateout[i,k,ll] = stateinit[i,k,ll] + dt * tend[i,k,ll]
            end
        end
    end
end

function set_halo_values_x!(state::OffsetArray{Float64, 3, Array{Float64, 3}})
        
end

function compute_tendencies_x!(state::OffsetArray{Float64, 3, Array{Float64, 3}},
                               flux::Array{Float64, 3},
                               tend::Array{Float64, 3})                                  
end

function set_halo_values_z!(state::OffsetArray{Float64, 3, Array{Float64, 3}})
    
end
        
function compute_tendencies_z!(state::OffsetArray{Float64, 3, Array{Float64, 3}},
                               flux::Array{Float64, 3},
                               tend::Array{Float64, 3})                           
end
            
function reductions(state::OffsetArray{Float64, 3, Array{Float64, 3}},
                    hy_dens_cell::OffsetVector{Float64, Vector{Float64}},
                    hy_dens_theta_cell::OffsetVector{Float64, Vector{Float64}})
    
    mass, te, r, u, w, th, p, t, ke, le = [zero(Float64) for _ in 1:10] 
    glob = Array{Float64}(undef, 2)
    
    for k in 1:NZ
        for i in 1:NX
            r  =   state[i,k,ID_DENS] + hy_dens_cell[k]             # Density
            u  =   state[i,k,ID_UMOM] / r                           # U-wind
            w  =   state[i,k,ID_WMOM] / r                           # W-wind
            th = ( state[i,k,ID_RHOT] + hy_dens_theta_cell[k] ) / r # Potential Temperature (theta)
            p  = C0*(r*th)^GAMMA      # Pressure
            t  = th / (P0/p)^(RD/CP)  # Temperature
            ke = r*(u*u+w*w)          # Kinetic Energy
            ie = r*CV*t               # Internal Energy
            mass = mass + r            *DX*DZ # Accumulate domain mass
            te   = te   + (ke + r*CV*t)*DX*DZ # Accumulate domain total energy
        end
    end
    
    MPI.Allreduce!(Array{Float64}([mass,te]), glob, +, COMM)
    
    return glob
end

function output()
    
end

function finalize!(state::OffsetArray{Float64, 3, Array{Float64, 3}})

    #println(axes(state))
    
end

# invoke main function
main(ARGS)
#@run main(ARGS)
