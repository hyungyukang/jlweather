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

JULIA := julia -g 2

SRCDIR := ../../src
JLSRC  := ${SRCDIR}/miniWeather_mpi.jl
ARGS := ${SIM_TIME} ${NX} ${NZ} ${DATA}

CXX := g++
FC  := gfortran

CXXFLAGS := -shared -fPIC

DATASRC := ${SRCDIR}/data.cpp
KERNELSRC := ${SRCDIR}/kernel.cpp

run:
	${JULIA} --color=yes -O -- ${JLSRC} ${ARGS}

slib: dlib.so klib.so

dlib.so: ${DATASRC}
	${CXX} ${CXXFLAGS} -o $@ $<
	mv $@ ${SRCDIR}

klib.so: ${KERNELSRC}
	${CXX} ${CXXFLAGS} -o $@ $<
	mv $@ ${SRCDIR}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

clean:
	rm -f *.so
