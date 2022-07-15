==================
Julia miniWeather
==================

This repo includes several versions of miniWeather:

==================== =============================================================================================================================
version              note 
==================== =============================================================================================================================
Fortran+MPI          Fortran MPI implementation copied from `original miniWeather repo <https://github.com/mrnorman/miniWeather/>`_
Fortran OpenAcc.     Fortran OpenAcc porting of miniWeather_mpi.F90 
Julia                Julia implementation of miniWeather_mpi.F90
Julia+JAI(OpenAcc)   Julia OpenAcc implementation using `JAI <https://github.com/grnydawn/AccelInterfaces.jl/>`_ of miniWeather_openacc.F90
Julia+Manual_OpenAcc Julia OpenAcc manual implementation of miniWeather_openacc.F90
==================== =============================================================================================================================


Fortran+MPI version
=====================

* src/fortran/miniWeather_mpi.F90

Fortran OpenAcc version
==========================

* src/fortran/miniWeather_openacc.F90

Julia version
==========================

* src/julia/miniWeather_mpi.jl

Julia+JAI(OpenAcc) version
==========================

* src/jai/miniWeather_accel.jl and .knl files

Julia+Manual_OpenAcc version
===============================

* src/manual/miniWeather_openacc.jl and Fortran driver files


..note::
	Running some of the versions create several output files such as shared libraries.
	For convinience, you may create a subdirectory under "run" directory for your experiment and run your command in the subdirectory.
	Please see **Makefile** in "run/crusher" subdirectory for more details.
