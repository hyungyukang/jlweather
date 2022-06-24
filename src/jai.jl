module AcceleratorInterface

#import Meta.parse

import Libdl.dlopen,
       Libdl.RTLD_LAZY,
       Libdl.RTLD_DEEPBIND,
       Libdl.RTLD_GLOBAL,
       Libdl.dlsym

import OffsetArrays.OffsetArray,
       OffsetArrays.OffsetVector

export AccelInfo, KernelInfo, copyin!, copyout!, launch!

struct AccelInfo
    
	sharedlib

	# load shared library
    function AccelInfo()

		run(`make datalib`)
		dlib = dlopen("./datalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
        new(dlib)
    end
end


struct KernelInfo
    
    ainfo::AccelInfo
	sharedlib
 
    function KernelInfo(ainfo::AccelInfo)

		run(`make kernellib`)
		klib = dlopen("./kernellib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
        new(ainfo, klib)
    end
end

# dynamically dispatch data to one of copyin functions generated in AccelInfo
# based on data...
function copyin!(ainfo::AccelInfo, data...)

    dataenter = dlsym(ainfo.sharedlib, :dataenter)

    args = []
    dtypes = []

    push!(args, length(data))
    push!(dtypes, typeof(args[end]))

    for arg in data
        #println(typeof(arg))
        if typeof(arg) <: OffsetArray

            offsets = arg.offsets

            push!(args, length(offsets))
            push!(dtypes, typeof(args[end]))

            push!(args, reverse(offsets))
            push!(dtypes, Ref{typeof(args[end])})

            push!(args, reverse(size(arg)))
            push!(dtypes, Ref{typeof(args[end])})

            push!(args, arg.parent)
            push!(dtypes, Ptr{typeof(args[end])})

            # for debugging
            #fill!(arg.parent, 100.)
            #arg.parent[end, end, end] = 100.
            arg.parent[3, 2, 1] = 100.

        elseif typeof(arg) <: AbstractArray

            offsets = Tuple(1 for x = 1:x)

            push!(args, length(offsets)) 
            push!(dtypes, typeof(args[end]))

            push!(args, reverse(offsets))
            push!(dtypes, Ref{typeof(args[end])})

            push!(args, reverse(size(arg)))
            push!(dtypes, Ref{typeof(args[end])})

            push!(args, arg)
            push!(dtypes, Ptr{typeof(args[end])})

            # for debugging
            #fill!(arg, 100.)

        else

            push!(args, Int64(0)) 
            push!(dtypes, typeof(args[end]))

            push!(args, arg)
            push!(dtypes, typeof(arg))
        end
    end

    argtypes = Meta.parse(string(((dtypes...),)))
    ccallexpr = :(ccall($dataenter, Int64, $argtypes, $(args...)))

    #val = ccall(dataenter, Int64, (Array{Float64, 3},), args[1])
    #@eval val = ccall($dataenter, Int64, (Array{Float64, 3},), $args[1])
    #@eval val = ccall($dataenter, Int64, (($dtypes...),), $args[1])
    #@eval val = ccall($dataenter, Int64, (($argtypes...),), (($args...),))
    #@eval val = ccall($dataenter, Int64, $argtypes, 3, ($aaa).offsets, size($aaa), ($aaa).parent )
    #@eval val = ccall($dataenter, Int64, $argtypes, $args[1], $args[2], $args[3], $args[4] )
    #@eval val = ccall($dataenter, Int64, $argtypes, $(argdata...))
    
    @eval val = $ccallexpr

    @show "CCC", val, data[1].parent[3,2,1]

    #if val == C_NULL
    #    error("dataenter: undefined variable: ", val)
    #end

end


# dynamically dispatch data to one of copyout functions generated in AccelInfo
# based on data...
function copyout!(ainfo::AccelInfo, data...)
end


# dynamically dispatch kernel one of launch functions generated in KernelInfo
# based on data...
function launch!(kinfo::KernelInfo, data...)
	println(kinfo)
end


# In case of multiprocessing, let master finish its work
function barrier(barrier_func, barrier_args...)
end

end
