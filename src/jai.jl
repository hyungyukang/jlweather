module AcceleratorInterface

using Debugger

import IniFile.Inifile

export AccelInfo, KernelInfo, copyto, copyfrom, launch

# static accelerator info to generate a shared library
# generate multiple copyto, copyfrom, and launch functions 
# based on union of data specifications
struct AccelInfo
    
    function AccelInfo()
        new()
    end
end


# static kernel info to generate a shared library
struct KernelInfo
    
    ainfo::AccelInfo
    kfile::AbstractString
    
    function KernelInfo(info::AccelInfo, inipath)

        ini = read(Inifile(), inipath)
        @bp
        
        new(info, inipath)
    end
end

# dynamically dispatch data to one of copyto functions generated in AccelInfo
# based on data...
function copyto(ainfo::AccelInfo, data...)

end


# dynamically dispatch data to one of copyfrom functions generated in AccelInfo
# based on data...
function copyfrom(ainfo::AccelInfo, data...)
end


# dynamically dispatch kernel one of launch functions generated in KernelInfo
# based on data...
function launch(kinfo::KernelInfo, data...)
end


# In case of multiprocessing, let master finish its work
function barrier(barrier_func, barrier_args...)
end

#dlib = Libdl.dlopen(PATH_DATALIB, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
#klib = Libdl.dlopen(PATH_KERNELLIB, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

#dataenter = Libdl.dlsym(dlib, :dataenter)
#val = ccall(dataenter, Int64, (Array{Float64, 3},), state.parent)
#@show "BBB", state.offsets

#if val == C_NULL
#    error("dataenter: undefined variable: ", val)
#end

end
