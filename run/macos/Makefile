SIM_TIME := 400.0
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

JULIA := julia --project=. -g 2

SRCDIR := ../../src
JLSRC  := ${SRCDIR}/miniWeather_mpi.jl
#ARGS := ${SIM_TIME} ${NX} ${NZ} ${DATA}
ARGS := 

CC := gcc
CXX := g++
FC  := gfortran

CXXFLAGS := -shared -fPIC
C_FLAGS := -fPIC -shared

run:
	${JULIA} --color=yes -O -- ${JLSRC} ${ARGS}


datalib: ${SRCDIR}/data.cpp
	${CC} ${C_FLAGS} -o $@.so $<

kernellib: ${SRCDIR}/kernel.cpp
	${CC} ${C_FLAGS} -o $@.so $<

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

clean:
	rm -f *.so
