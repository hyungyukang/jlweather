Julia miniWeather
==================

This repo includes several versions of miniWeather:

* src/fortran/miniWeather_mpi.F90: Fortran MPI implementation copied from `original miniWeather repo <https://github.com/mrnorman/miniWeather/>`_
* src/fortran/miniWeather_openacc.F90: Fortran OpenAcc porting of miniWeather_mpi.F90 
* src/julia/miniWeather_mpi.jl: Julia implementation of miniWeather_mpi.F90
* src/jai/miniWeather_accel.jl and .knl files: Julia OpenAcc implementation using `JAI <https://github.com/grnydawn/AccelInterfaces.jl/>`_ of miniWeather_openacc.F90
* src/manual/miniWeather_openacc.jl and driver fortran files: Julia OpenAcc manual implementation of miniWeather_openacc.F90

Running the version creates several files such as shared libraries. As a convention, you may create a subdirectory under "run" directory for your experiment and run your command in the directory. Please see Makefiles in pre-existing the subdirectoris for more details.
