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
    
    acceltype
	sharedlib

	# load shared library
    function AccelInfo(atype::String="F")

        if atype == "C"
            run(`make c_datalib`)
            dlib = dlopen("./c_datalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        elseif atype == "F"
            run(`make f_datalib`)
            dlib = dlopen("./f_datalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        end

        new(atype, dlib)
    end
end


struct KernelInfo
    
    ainfo::AccelInfo
	sharedlib
 
    function KernelInfo(ainfo::AccelInfo)

        if ainfo.acceltype == "C"

		    run(`make c_kernellib`)
		    klib = dlopen("./c_kernellib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        elseif ainfo.acceltype == "F"

		    run(`make f_kernellib`)
		    klib = dlopen("./f_kernellib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        end

        new(ainfo, klib)
    end
end

# dynamically dispatch data to one of copyin functions generated in AccelInfo
# based on data...
function copyin!(ainfo::AccelInfo, data...)

    dataenter = dlsym(ainfo.sharedlib, :dataenter)

    args = []
    dtypes = []

    if ainfo.acceltype == "C"

        push!(args, length(data))
        push!(dtypes, typeof(args[end]))

    elseif ainfo.acceltype == "F"

        push!(args, length(data))
        push!(dtypes, Ref{typeof(args[end])})

    end


    for arg in data
        #println(typeof(arg))
        
        if typeof(arg) <: OffsetArray

            offsets = arg.offsets

            if ainfo.acceltype == "C"
                push!(args, length(offsets))
                push!(dtypes, typeof(args[end]))

                push!(args, reverse(offsets))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, reverse(size(arg)))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg.parent)
                push!(dtypes, Ptr{typeof(args[end])})

                # for debugging
                #arg.parent[end, end, end] = 100.
                arg.parent[3, 2, 1] = 100.

            elseif ainfo.acceltype == "F"
                push!(args, length(offsets))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, offsets)
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, size(arg))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg.parent)
                push!(dtypes, Ptr{typeof(args[end])})

                # for debugging
                arg.parent[3, 2, 1] = 100.

            end

        elseif typeof(arg) <: AbstractArray

            offsets = Tuple(1 for x = 1:x)

            if ainfo.acceltype == "C"

                push!(args, length(offsets)) 
                push!(dtypes, typeof(args[end]))

                push!(args, reverse(offsets))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, reverse(size(arg)))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg)
                push!(dtypes, Ptr{typeof(args[end])})

            elseif ainfo.acceltype == "F"

                push!(args, length(offsets)) 
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, offsets)
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, size(arg))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg)
                push!(dtypes, Ptr{typeof(args[end])})

            end
        else

            if ainfo.acceltype == "C"

                push!(args, Int64(0)) 
                push!(dtypes, typeof(args[end]))

                push!(args, arg)
                push!(dtypes, typeof(args[end]))

            elseif ainfo.acceltype == "F"

                push!(args, Int64(0)) 
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg)
                push!(dtypes, Ref{typeof(args[end])})

            end
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
