module AcceleratorInterface

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

	dtypes = ()
	args = []

    for arg in data
        #println(typeof(arg))
        if arg isa OffsetArray

            println(typeof(arg.parent), typeof(arg.offsets))
            push!(dtypes, typeof(arg.parent))
            push!(args, arg.parent)


        elseif arg isa AbstractArray
            push!(dtypes, typeof(arg))
            push!(args, arg)
        else
            push!(dtypes, typeof(arg))
            push!(args, arg)
        end
    end

    println("BBBB", typeof(dtypes), typeof(args))
    #val = ccall(dataenter, Int64, dtypes, 1)

    #@show "BBB", state.offsets

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
