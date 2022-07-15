
Julia miniWeather
==================

This repo includes several versions of miniWeather.

+----------------------+------------------------------------------------------------------------------------------------------------------------+
|      version         |          note                                                                                                          |
+======================+========================================================================================================================+
| Fortran MPI          | Fortran MPI implementation copied from `original miniWeather repo <https://github.com/mrnorman/miniWeather/>`_         |
+----------------------+------------------------------------------------------------------------------------------------------------------------+
| Fortran OpenAcc      | Fortran OpenAcc porting of Fortran MPI version                                                                         |
+----------------------+------------------------------------------------------------------------------------------------------------------------+
| Julia                | Julia implementation of Fortran MPI version                                                                            |
+----------------------+------------------------------------------------------------------------------------------------------------------------+
| Julia Manual OpenAcc | Julia OpenAcc manual implementation of Fortran OpenAcc version                                                         |
+----------------------+------------------------------------------------------------------------------------------------------------------------+
| Julia JAI(OpenAcc)   | Julia OpenAcc implementation using `JAI <https://github.com/grnydawn/AccelInterfaces.jl/>`_ of Fortran OpenAcc version |
+----------------------+------------------------------------------------------------------------------------------------------------------------+

Fortran MPI version
=====================

* src/fortran/miniWeather_mpi.F90

Fortran OpenAcc version
==========================

* src/fortran/miniWeather_openacc.F90

Julia version
==========================

* src/julia/miniWeather_mpi.jl

Julia Manual OpenAcc version
===============================

* src/manual/miniWeather_openacc.jl and Fortran driver files

Julia JAI(OpenAcc) version
==========================

* src/jai/miniWeather_accel.jl and .knl files

**Notes**
    Running some of the versions create several output files such as shared libraries.
    For convinience, you may create a subdirectory under "run" directory for your experiment and run your command in the subdirectory.
    Please see **Makefile** in "run/crusher" subdirectory for more details.
