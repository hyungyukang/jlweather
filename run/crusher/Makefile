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

#DATA := ${DATA_SPEC_COLLISION}
DATA := ${DATA_SPEC_THERMAL}

#JULIA := julia --project=.
JULIA := julia --project=/ccs/home/grnydawn/repos/github/jlweather/run/crusher
SRCDIR :=../../src
SRC := ${SRCDIR}/miniWeather_mpi.jl
ARGS := -s ${SIM_TIME} -x ${NX} -z ${NZ} -f ${OUT_FREQ} -d ${DATA}
#ARGS := 

CC := cc
CXX := CC
FC := ftn

CXX_FLAGS := -fPIC -shared
C_FLAGS := -fPIC -shared
F_FLAGS := -fPIC -shared

run:
	${JULIA} --color=yes -O -- ${SRC} ${ARGS}

c_enterdatalib: ${SRCDIR}/enterdata.cpp
	${CC} ${C_FLAGS} -o $@.so $<

c_exitdatalib: ${SRCDIR}/exitdata.cpp
	${CC} ${C_FLAGS} -o $@.so $<

c_kernellib: ${SRCDIR}/kernel.cpp
	${CC} ${C_FLAGS} -o $@.so $<

f_enterdatalib: ${SRCDIR}/enterdata.F90
	${FC} ${F_FLAGS} -o $@.so $<

f_exitdatalib: ${SRCDIR}/exitdata.F90
	${FC} ${F_FLAGS} -o $@.so $<

f_kernellib: ${SRCDIR}/kernel.F90
	${FC} ${F_FLAGS} -o $@.so $<

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

clean:
	rm -rf *.so core *.mod .jaitmp *.o *.s
