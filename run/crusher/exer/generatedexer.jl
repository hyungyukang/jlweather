# NOTES
###########

# question 1: what is the best way to generate function according to the type of arguments AND the size of array arguments

# questions2: if there are many generated fuctions withe the same function name, say, several thousants, how much impact to performance

# thoughts: I may generate different function names for differnt type and size, and compiler... of array argument

int1 = Vector{Int64}(undef, 1)
int2 = Array{Int64, 2}(undef, 2, 3)
float1 = Vector{Float64}(undef, 3)
float2 = Array{Float64, 2}(undef, 2, 3)

println(int1, int2, float1, float2)

@generated function ex1(expr...)
    println("POS1 : ", expr...)
    #return :(println("POS2 : ", $(expr...)); expr...)
    return :(println("POS2 : ", $(expr...)); expr)
end


println("OUT1 :", ex1(int1))
println("OUT2 :", ex1(int2))
println("OUT3 :", ex1(int1, int2))

println("OUT4 :", ex1(float1))
println("OUT5 :", ex1(float2))
println("OUT6 :", ex1(float1, float2))
