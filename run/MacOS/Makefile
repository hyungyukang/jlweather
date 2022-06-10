SIM_TIME := 1000.0
NX := 2000
NZ := 100

DATA_SPEC_COLLISION       := 1
DATA_SPEC_THERMAL         := 2
DATA_SPEC_MOUNTAIN        := 3
DATA_SPEC_TURBULENCE      := 4
DATA_SPEC_DENSITY_CURRENT := 5
DATA_SPEC_INJECTION       := 6

DATA := ${DATA_SPEC_COLLISION}

JULIA := julia
SRC := ../../src/miniWeather_mpi.jl
ARGS := ${SIM_TIME} ${NX} ${NZ} ${DATA}

run:
	${JULIA} --color=yes -O -- ${SRC} ${ARGS}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 ${SRC} ${ARGS}

