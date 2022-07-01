function ex1(test)
    println("in ex1: test = ", test)
end

ex1(1)

@generated function ex2(expr1)
    println("in compile ex2: expr1 = ", expr1)
    return :(println("in execute ex2: ", $expr1); expr1)
end


println("WWW", ex2(1+1))
println("WWW", ex2(1+1))
println("WWW", ex2(1.0+1.0))
println("WWW", ex2("test"))
println("WWW", ex2("test"))
println("WWW", ex2("test1"))
