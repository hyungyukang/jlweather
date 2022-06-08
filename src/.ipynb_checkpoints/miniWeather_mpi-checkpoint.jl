using OffsetArrays

##############
# constants
##############

const SIM_TIME    = 1000.0
const NX_GLOB     = 100
const NZ_GLOB     = 50
const HS          = 4
const NUM_VARS    = 4
const XLEN        = Float64(2.e4) # Length of the domain in the x-direction (meters)
const ZLEN        = Float64(1.e4) # Length of the domain in the z-direction (meters)
const CFL         = Float64(1.5)  # "Courant, Friedrichs, Lewy" number (for numerical stability)
const MAX_SPEED   = Float64(450.0)# Assumed maximum wave speed during the simulation (speed of sound + speed of wind) (meter / sec)
const DX          = XLEN / NX_GLOB
const DZ          = ZLEN / NZ_GLOB
const DT          = min(DX,DZ) / MAX_SPEED * CFL

        
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
function main()

    ######################
    # top-level variables
    ######################

    etime = Float64(0.0)

    dt, state, statetmp, flux, tend = init!()

    reductions()

    output()

    # main loop
    while etime < SIM_TIME

        if etime + dt > SIM_TIME
            dt = SIM_TIME - etime
        end

        timestep(state, statetmp, flux, tend, dt)

        etime = etime + dt

    end

    reductions()

    finalize!(state)

end

function init!()

    nranks = 1
    myrank = 0
    
    nper = Float64(NX_GLOB)/nranks
    i_beg = trunc(Int, round( nper* (myrank)    )+1)
    i_end = trunc(Int, round( nper*((myrank)+1) ))
    nx = i_end - i_beg + 1
    left_rank  = myrank - 1
    if left_rank == -1
        left_rank = nranks-1
    end
    right_rank = myrank + 1
    if right_rank == nranks
        right_rank = 0
    end
            
    nz = NZ_GLOB
    
    #dt, state, statetmp, flux, tend = init()
    dt = Float64(1.0)
    
    #_state      = Array{Float64}(undef, nx+2*HS, nz+2*HS, NUM_VARS) 
    _state      = zeros(Float64, nx+2*HS, nz+2*HS, NUM_VARS) 
    state       = OffsetArray(_state, 1-HS:nx+HS, 1-HS:nz+HS, 1:NUM_VARS)
    _statetmp   = Array{Float64}(undef, nx+2*HS, nz+2*HS, NUM_VARS) 
    statetmp    = OffsetArray(_state, 1-HS:nx+HS, 1-HS:nz+HS, 1:NUM_VARS)
    flux        = Array{Float64}(undef, nx+1, nz+1, NUM_VARS) 
    tend        = Array{Float64}(undef, nx, nz, NUM_VARS) 

    return dt, state, statetmp, flux, tend
end

function timestep(state, statetmp, flux, tend, dt)
    
end

        
function reductions()
    
end

function output()
    
end

function finalize!(state)

    #println(axes(state))
    
end

# invoke main function
main()
