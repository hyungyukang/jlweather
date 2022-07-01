module AccelInterfaces

include("./flang.jl")
import  .FortranInterface as JFI

import Libdl.dlopen,
       Libdl.RTLD_LAZY,
       Libdl.RTLD_DEEPBIND,
       Libdl.RTLD_GLOBAL,
       Libdl.dlsym

import OffsetArrays.OffsetArray,
       OffsetArrays.OffsetVector

export AccelInfo, KernelInfo, copyin!, copyout!, launch!,
       FLANG, CLANG

@enum AccelType FLANG CLANG

struct AccelInfo
    
    acceltype::AccelType
    ismaster::Bool
    sharedlibs::Dict

    function AccelInfo(acceltype::AccelType=FLANG; ismaster::Bool=true)

        new(acceltype, ismaster, Dict())
    end
end

struct KernelInfo
    
    kernelpath::String

    function KernelInfo(path::String)

        new(path)
    end
end


function argsdtypes(ainfo::AccelInfo, copyin, copyout, copyinout)

    args = []
    dtypes = []
    data = []

    for cin in copyin
        push!(data, cin)
    end
    for cout in copyout
        push!(data, cout)
    end
    for cinout in copyinout
        push!(data, cinout)
    end

    for arg in data
        
        if typeof(arg) <: OffsetArray

            offsets = arg.offsets

            if ainfo.acceltype == CLANG

                push!(args, arg.parent)
                push!(dtypes, Ptr{typeof(args[end])})

            elseif ainfo.acceltype == FLANG

                # for debugging
                #arg[0, 0, 1] = 100.
                #arg.parent[2, 2, 1] = 100.

                #push!(args, length(offsets))
                #push!(dtypes, Ref{typeof(args[end])})

                #push!(args, offsets)
                #push!(dtypes, Ref{typeof(args[end])})

                #push!(args, size(arg))
                #push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg.parent)
                push!(dtypes, Ptr{typeof(args[end])})

            end

        elseif typeof(arg) <: AbstractArray

            #offsets = Tuple(1 for x = 1:x)

            #push!(args, length(offsets)) 
            #push!(dtypes, typeof(args[end]))

            #push!(args, reverse(offsets))
            #push!(dtypes, Ref{typeof(args[end])})

            #push!(args, reverse(size(arg)))
            #push!(dtypes, Ref{typeof(args[end])})

            push!(args, arg)
            push!(dtypes, Ptr{typeof(args[end])})

        else

            if ainfo.acceltype == CLANG

                push!(args, arg)
                push!(dtypes, typeof(args[end]))

            elseif ainfo.acceltype == FLANG

                push!(args, arg)
                push!(dtypes, Ref{typeof(args[end])})

            end
        end
    end

    args, dtypes
end

function argsdtypes1(ainfo::AccelInfo, data...)

    args = []
    dtypes = []

    if ainfo.acceltype == CLANG

        push!(args, length(data))
        push!(dtypes, typeof(args[end]))

    elseif ainfo.acceltype == FLANG

        push!(args, length(data))
        push!(dtypes, Ref{typeof(args[end])})

    end

    for arg in data
        #println(typeof(arg))
        
        if typeof(arg) <: OffsetArray

            offsets = arg.offsets

            if ainfo.acceltype == CLANG
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

            elseif ainfo.acceltype == FLANG
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

            if ainfo.acceltype == CLANG

                push!(args, length(offsets)) 
                push!(dtypes, typeof(args[end]))

                push!(args, reverse(offsets))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, reverse(size(arg)))
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg)
                push!(dtypes, Ptr{typeof(args[end])})

            elseif ainfo.acceltype == FLANG

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

            if ainfo.acceltype == CLANG

                push!(args, Int64(0)) 
                push!(dtypes, typeof(args[end]))

                push!(args, arg)
                push!(dtypes, typeof(args[end]))

            elseif ainfo.acceltype == FLANG

                push!(args, Int64(0)) 
                push!(dtypes, Ref{typeof(args[end])})

                push!(args, arg)
                push!(dtypes, Ref{typeof(args[end])})

            end
        end
    end

    args, dtypes
end


function copyin!(ainfo::AccelInfo, data...)

    # generate sha from data

    args, dtypes = argsdtypes(ainfo, data, (), ())

    # generate hash for enterdata

    datahash = hash(("copyin!", ainfo.acceltype, dtypes))

    if !haskey(ainfo.sharedlibs, datahash)
        println("Compiling shared library")

        if ainfo.acceltype == CLANG
            run(`make c_enterdatalib`)
            dlib = dlopen("./c_enterdatalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        elseif ainfo.acceltype == FLANG

            # generate source code
            #
            # compile code

            run(`make f_enterdatalib`)
            dlib = dlopen("./f_enterdatalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        end

        ainfo.sharedlibs[datahash] = dlib

    else
        println("Reusing shared library")
        dlib = ainfo.sharedlibs[datahash]

    end

    dataenter = dlsym(dlib, :dataenter)

    argtypes = Meta.parse(string(((dtypes...),)))
    # C and Fortran does not need to enter data
    #ccallexpr = :(ccall($dataenter, Int64, $argtypes, $(args...)))
    ccallexpr = :(ccall($dataenter, Int64, ()))

    #val = ccall(dataenter, Int64, (Array{Float64, 3},), args[1])
    #@eval val = ccall($dataenter, Int64, (Array{Float64, 3},), $args[1])
    #@eval val = ccall($dataenter, Int64, (($dtypes...),), $args[1])
    #@eval val = ccall($dataenter, Int64, (($argtypes...),), (($args...),))
    #@eval val = ccall($dataenter, Int64, $argtypes, 3, ($aaa).offsets, size($aaa), ($aaa).parent )
    #@eval val = ccall($dataenter, Int64, $argtypes, $args[1], $args[2], $args[3], $args[4] )
    #@eval val = ccall($dataenter, Int64, $argtypes, $(argdata...))
    
    @eval val = $ccallexpr

    #@show "CCC", val, data[1].parent[3,2,1]

    #if val == C_NULL
    #    error("dataenter: undefined variable: ", val)
    #end

end


# dynamically dispatch data to one of copyout functions generated in AccelInfo
# based on data...
function copyout!(ainfo::AccelInfo, data...)

    # generate sha from data

    args, dtypes = argsdtypes(ainfo, (), data, ())

    # generate hash for exitdata

    datahash = hash(("copyout!", ainfo.acceltype, dtypes))

    if !haskey(ainfo.sharedlibs, datahash)
        println("Compiling shared library")

        if ainfo.acceltype == CLANG
            run(`make c_exitdatalib`)
            dlib = dlopen("./c_exitdatalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        elseif ainfo.acceltype == FLANG

            # generate source code
            #
            # compile code

            run(`make f_exitdatalib`)
            dlib = dlopen("./f_exitdatalib.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        end

        ainfo.sharedlibs[datahash] = dlib

    else
        println("Reusing shared library")
        dlib = ainfo.sharedlibs[datahash]

    end

    dataexit = dlsym(dlib, :dataexit)

    argtypes = Meta.parse(string(((dtypes...),)))
    # C and Fortran does not need to exit data
    #ccallexpr = :(ccall($dataexit, Int64, $argtypes, $(args...)))
    ccallexpr = :(ccall($dataexit, Int64, ()))
    
    @eval val = $ccallexpr

    @show "CCC", val, data[1].parent[3,2,1]

end


# dynamically dispatch kernel one of launch functions generated in 
# based on data...
function launch!(ainfo::AccelInfo, kinfo::KernelInfo,
                data...; copyout::Union{Tuple, Vector}=(),
                copyinout::Union{Tuple, Vector}=(),
                compile::Union{String, Nothing}=nothing,
                workdir::Any=nothing)

    # generate sha from data

    args, dtypes = argsdtypes(ainfo, data, copyout, copyinout)

    # generate hash for exitdata

    kernelhash = hash(("launch!", ainfo.acceltype, dtypes))

    if !haskey(ainfo.sharedlibs, kernelhash)
        println("Compiling shared library")

        if ainfo.acceltype == CLANG
            run(`make c_$kernelhash KERNELPATH=$kinfo.kernelpath`)
            dlib = dlopen("./c_$kernelhash.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        elseif ainfo.acceltype == FLANG

            # generate source code
            #
            # compile code
            
            fname, ext = splitext(basename(kinfo.kernelpath))

            if compile === nothing
                run(`make f_$fname KERNELPATH=$(kinfo.kernelpath)`)

            else
                run(`$(split(compile)) -o f_$fname.so $(kinfo.kernelpath)`)

            end

            dlib = dlopen("./f_$fname.so", RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)

        end

        ainfo.sharedlibs[kernelhash] = dlib

    else
        println("Reusing shared library")
        dlib = ainfo.sharedlibs[kernelhash]

    end

    launch = dlsym(dlib, :launch)

    argtypes = Meta.parse(string(((dtypes...),)))
    # C and Fortran does not need to exit data
    ccallexpr = :(ccall($launch, Int64, $argtypes, $(args...)))
    
    @eval val = $ccallexpr

end


# In case of multiprocessing, let master finish its work
function barrier(barrier_func, barrier_args...)
end

end
