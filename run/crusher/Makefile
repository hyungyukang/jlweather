# example command: jsrun  -n 2 -- make
# example command: srun  -n 2 -- make

SIM_TIME := 400.0
OUT_FREQ := 100.0
NX := 100 
NZ := 50
#SIM_TIME := 10.0
#OUT_FREQ := 10.0
#NX := 640 
#NZ := 100

DATA_SPEC_COLLISION       := 1
DATA_SPEC_THERMAL         := 2
DATA_SPEC_MOUNTAIN        := 3
DATA_SPEC_TURBULENCE      := 4
DATA_SPEC_DENSITY_CURRENT := 5
DATA_SPEC_INJECTION       := 6

#DATA_SPEC := ${DATA_SPEC_COLLISION}
DATA_SPEC := ${DATA_SPEC_THERMAL}

#JULIA := julia --project=.
JULIA := julia --project=/ccs/home/grnydawn/repos/github/jlweather/run/crusher
SRCDIR :=../../src
SRC := ${SRCDIR}/miniWeather_mpi.jl
FORTSRC := ${SRCDIR}/miniWeather_mpi.F90
ARGS := -s ${SIM_TIME} -x ${NX} -z ${NZ} -f ${OUT_FREQ} -d ${DATA_SPEC}
#ARGS := 

CC := cc
CXX := CC
FC := ftn

INCLUDES := -I${OLCF_PARALLEL_NETCDF_ROOT}/include
LIBS := -L${OLCF_PARALLEL_NETCDF_ROOT}/lib -lpnetcdf
MACROS := -D_NX=${NX} -D_NZ=${NZ} -D_SIM_TIME=${SIM_TIME} -D_OUT_FREQ=${OUT_FREQ} -D_DATA_SPEC=${DATA_SPEC}

F_FLAGS := ${INCLUDES} ${LIBS} ${MACROS}

julia:
	${JULIA} --color=yes -O -- ${SRC} ${ARGS}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

fort: ${FORTSRC}
	${FC} ${F_FLAGS} -o miniweather_fort.exe $<

clean:
	rm -rf *.so core *.mod .jaitmp *.o *.s miniweather_fort.exe output.nc
