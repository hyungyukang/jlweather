# example command: jsrun  -n 2 -- make

SIM_TIME := 400.0
OUT_FREQ := 100.0
NX := 100 
NZ := 50

DATA_SPEC_COLLISION       := 1
DATA_SPEC_THERMAL         := 2
DATA_SPEC_MOUNTAIN        := 3
DATA_SPEC_TURBULENCE      := 4
DATA_SPEC_DENSITY_CURRENT := 5
DATA_SPEC_INJECTION       := 6

#DATA := ${DATA_SPEC_COLLISION}
DATA := ${DATA_SPEC_THERMAL}

JULIA := julia --project=.
SRCDIR :=../../src
SRC := ${SRCDIR}/miniWeather_mpi.jl
#ARGS := ${SIM_TIME} ${NX} ${NZ} ${OUT_FREQ} ${DATA}
ARGS := 
C_FLAGS := -fPIC -shared
CC := cc

run:
	${JULIA} --color=yes -O -- ${SRC} ${ARGS}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

datalib: ${SRCDIR}/data.cpp
	${CC} ${C_FLAGS} -o $@.so $<

kernellib: ${SRCDIR}/kernel.cpp
	${CC} ${C_FLAGS} -o $@.so $<

clean:
	rm -rf core output.nc *.so
