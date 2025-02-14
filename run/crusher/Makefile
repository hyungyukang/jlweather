# example command: srun  -n 2 -- make

#SIM_TIME := 400.0
OUT_FREQ := 100.0
SIM_TIME := 100.0
NX := 200
NZ := 100
#SIM_TIME := 10.0 #10.0
#OUT_FREQ := 10.0 # 10.0
#NX := 6400 # 640 
#NZ := 100 # 100

DATA_SPEC_COLLISION       := 1
DATA_SPEC_THERMAL         := 2
DATA_SPEC_MOUNTAIN        := 3
DATA_SPEC_TURBULENCE      := 4
DATA_SPEC_DENSITY_CURRENT := 5
DATA_SPEC_INJECTION       := 6

#DATA_SPEC := ${DATA_SPEC_COLLISION}
DATA_SPEC := ${DATA_SPEC_THERMAL}

OUT_FILE := /gpfs/alpine/cli133/proj-shared/grnydawn/temp/output.nc

JULIA := julia --project=.

SRCDIR :=../../src
JAISRC := ${SRCDIR}/jai/miniWeather_accel.jl
JLSRC := ${SRCDIR}/julia/miniWeather_mpi.jl
FORTSRC := ${SRCDIR}/fortran/miniWeather_mpi.F90
FORTACCSRC := ${SRCDIR}/fortran/miniWeather_openacc.F90
MANUALJLSRC := ${SRCDIR}/manual/miniWeather_openacc.jl
MANUALFORTSRC := ${SRCDIR}/manual/openacc_driver_${NX}x${NZ}.F90

ARGS := -s ${SIM_TIME} -x ${NX} -z ${NZ} -f ${OUT_FREQ} -d ${DATA_SPEC} -o ${OUT_FILE}

CC := cc
CXX := CC
FC := ftn

INCLUDES := -I${OLCF_PARALLEL_NETCDF_ROOT}/include
LIBS := -L${OLCF_PARALLEL_NETCDF_ROOT}/lib -lpnetcdf
MACROS := -D_NX=${NX} -D_NZ=${NZ} -D_SIM_TIME=${SIM_TIME} -D_OUT_FREQ=${OUT_FREQ} -D_DATA_SPEC=${DATA_SPEC}

F_FLAGS := ${INCLUDES} ${LIBS} ${MACROS}

jai:
	${JULIA} --color=yes -O -- ${JAISRC} ${ARGS}

julia:
	${JULIA} --color=yes -O -- ${JLSRC} ${ARGS}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

fort: ${FORTSRC}
	${FC} ${F_FLAGS} -o miniweather_fort.exe $<

fortacc: ${FORTACCSRC}
	${FC} ${F_FLAGS} -h acc,noomp -o miniweather_fortacc.exe $<

manual: ${MANUALFORTSRC}
	ftn -shared -fPIC -h acc,noomp -o jopenacc.so $<
	${JULIA} ${MANUALJLSRC} ${ARGS}

clean:
	rm -rf *.so core *.mod .jaitmp *.o *.s *.exe output.nc
