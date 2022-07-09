macro split(funcs...)
    tmp = Expr(:block)
    for func in funcs 
        kwexpr = Expr(:kw, :src, __source__)
        push!(func.args, kwexpr)
        dump(func)
        push!(tmp.args, func)    
    end
    return(tmp)
end


function f1(data; src::Any=nothing)
    println("This is f1.", data, src)
end

function f2(data; src::Any=nothing)
    println("This is f2.", data, src)
end

@split f1(1) f2(2)
