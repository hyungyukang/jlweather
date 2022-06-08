
JULIA := julia
SRC := ../../src/miniWeather_mpi.jl

run:
	${JULIA} --color=yes -O -- ${SRC}

debug:
	${JULIA} --color=yes --check-bounds=yes -O0 -- ${SRC}

