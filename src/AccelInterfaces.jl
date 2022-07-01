module AccelInterfaces

export AccelInfo, KernelInfo, FLANG, CLANG, get_accel!,get_kernel!,
       allocate!, deallocate!, copyin!, copyout!, launch!

@enum AccelType FLANG CLANG ANYACCEL

struct AccelInfo
    
    acceltype::AccelType
    ismaster::Bool
    constvars::Tuple
    constnames::NTuple
    compile::Union{String, Nothing}
    sharedlibs::Dict
    constants::Dict

    function AccelInfo(acceltype::AccelType; ismaster::Bool=true, constvars::Tuple,
                    compile::Union{String, Nothing}=nothing,
                    constnames::NTuple=())

        new(acceltype, ismaster, constvars, constnames, compile, Dict(), Dict())
    end
end


struct KernelInfo
    
    accel::AccelInfo
    kernelpath::String

    function KernelInfo(accel::AccelInfo, path::String)

        new(accel, path)
    end
end

function get_accel!(acceltype::AccelType; ismaster::Bool=true, constvars::Tuple=(),
                    compile::Union{String, Nothing}=nothing,
                    constnames::NTuple=())

    return AccelInfo(acceltype, ismaster=ismaster, constvars=constvars,
                    compile=compile, constnames=constnames)
end

function get_kernel!(accel::AccelInfo, path::String)

    return KernelInfo(accel, path)
end

function allocate!(accel::AccelInfo, data...)
end

function allocate!(kernel::KernelInfo, data...)
    return allocate!(kernel.accel, data...)
end

function deallocate!(accel::AccelInfo, data...)
end

function deallocate!(kernel::KernelInfo, data...)
    return deallocate!(kernel.accel, data...)
end

function copyin!(accel::AccelInfo, data...)
end

function copyin!(kernel::KernelInfo, data...)
    return copyin!(kernel.accel, data...)
end

function copyout!(accel::AccelInfo, data...)
end

function copyout!(kernel::KernelInfo, data...)
    return copyout!(kernel.accel, data...)
end

function launch!(kernel::KernelInfo, invars...;
                 innames::NTuple=(), outnames=NTuple=(),
                 outvars::Union{Tuple, Vector}=(),
                 compile::Union{String, Nothing}=nothing)

	println("TTT", length(kernel.accel.constnames))
end

end
