module AcceleratorInterface

include("./kernel.jl")
import  .Kernel.parse_kernel

#using Debugger

export AccelInfo, KernelInfo, copyin!, copyout!, launch!

# static accelerator info to generate a shared library
# generate multiple copyin, copyout, and launch functions 
# based on union of data specifications
struct AccelInfo
    
    function AccelInfo()
        new()
    end
end


# static kernel info to generate a shared library
struct KernelInfo
    
    ainfo::AccelInfo
    path::AbstractString
	parsed::Vector
 
    function KernelInfo(ainfo::AccelInfo, inipath, argtypes:Any)

		parsed = parse_kernel(inipath)

		# select accelerator from ainfo
		acceltype = "fortran"

		# generate code
		section = get_section(acceltype)

		# compile code
		libpath = Fortran.gen_library(section, argtypes)

		# load shared library

		#dlib = Libdl.dlopen(PATH_DATALIB, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
		klib = Libdl.dlopen(libpath, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

		#dataenter = Libdl.dlsym(dlib, :dataenter)
		#val = ccall(dataenter, Int64, (Array{Float64, 3},), state.parent)
		#@show "BBB", state.offsets

		#if val == C_NULL
		#    error("dataenter: undefined variable: ", val)
		#end

        new(ainfo, inipath, parsed)
    end
end

# dynamically dispatch data to one of copyin functions generated in AccelInfo
# based on data...
function copyin!(ainfo::AccelInfo, data...)

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
