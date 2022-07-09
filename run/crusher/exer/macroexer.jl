macro twostep(arg)
           println("I execute at parse time. The argument is: ", arg)
            println("__source__.file = ", __source__.file)
            println("__source__.line = ", __source__.line)
           return :(println("I execute at runtime. The argument is: ", $arg))
       end

ex = macroexpand(Main, :(@twostep :(1, 2, 3)) )

println("typeof(ex) = ", typeof(ex))

println("ex = ", ex)

println("eval(ex) = ", eval(ex))

@twostep(1)

