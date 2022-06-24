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

    dtypes = []
	args = []

    for arg in data
        #println(typeof(arg))
        if typeof(arg) <: OffsetArray

            offsets = arg.offsets

            push!(args, length(offsets))
            push!(dtypes, typeof(args[end]))

            push!(args, offsets)
            push!(dtypes, typeof(args[end]))

            push!(args, size(arg))
            push!(dtypes, typeof(args[end]))

            push!(args, arg.parent)
            push!(dtypes, typeof(args[end]))

        elseif typeof(arg) <: AbstractArray

            offsets = Tuple(1 for x = 1:x)

            push!(args, length(offsets)) 
            push!(dtypes, typeof(args[end]))

            push!(args, offsets)
            push!(dtypes, typeof(args[end]))

            push!(args, size(arg))
            push!(dtypes, typeof(args[end]))

            push!(args, arg)
            push!(dtypes, typeof(args[end]))

        else

            push!(args, Int64(0)) 
            push!(dtypes, typeof(args[end]))

            push!(args, arg)
            push!(dtypes, typeof(arg))
        end
    end

    strtypes = string(((dtypes...),))
    argtypes = Meta.parse(strtypes)

    strargs = "(args[1], args[2], args[3], args[4])"
    #println("BBB", strargs)
    argdata = Meta.parse(strargs)

    xx = :(ccall($dataenter, Int64, $argtypes, $(args...)))
    #println("XXX", xx)

    #println("CCC", dtypes, typeof(args))
    #val = ccall(dataenter, Int64, (Array{Float64, 3},), args[1])
    #@eval val = ccall($dataenter, Int64, (Array{Float64, 3},), $args[1])
    #@eval val = ccall($dataenter, Int64, (($dtypes...),), $args[1])
    #@eval val = ccall($dataenter, Int64, (($argtypes...),), (($args...),))
    #@eval val = ccall($dataenter, Int64, $argtypes, 3, ($aaa).offsets, size($aaa), ($aaa).parent )
    #@eval val = ccall($dataenter, Int64, $argtypes, $args[1], $args[2], $args[3], $args[4] )
    #@eval val = ccall($dataenter, Int64, $argtypes, $(argdata...))
    @eval val = $xx
    #val = ccall(dataenter, Int64, (Array{Float64, 3},), args...)

    @show "CCC", val

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
