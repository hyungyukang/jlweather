Julia miniWeather
==================

This repo includes several versions of miniWeather:

============================================= ======================================================================================================================================================================================
version                                        note 
============================================= ======================================================================================================================================================================================
src/fortran/miniWeather_mpi.F90.              Fortran MPI implementation copied from `original miniWeather repo <https://github.com/mrnorman/miniWeather/>`_
src/fortran/miniWeather_openacc.F90.          Fortran OpenAcc porting of miniWeather_mpi.F90 
src/julia/miniWeather_mpi.jl                  Julia implementation of miniWeather_mpi.F90
src/jai/miniWeather_accel.jl and .knl files   Julia OpenAcc implementation using `JAI <https://github.com/grnydawn/AccelInterfaces.jl/>`_ of miniWeather_openacc.F90
src/manual/miniWeather_openacc.jl and drivers Julia OpenAcc manual implementation of miniWeather_openacc.F90
============================================= ======================================================================================================================================================================================

Running some of the versions create several output files such as shared libraries. For convinience, you may create a subdirectory under "run" directory for your experiment and run your command in the subdirectory. Please see **Makefile** in "run/crusher" subdirectory for more details.
