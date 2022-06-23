module Kernel

# Kernel has sections
# Sections has Section

export parse_kernel

struct Section
    header::AbstractString
    body::Vector{AbstractString}
end

function parse_kernel(kernelfile::AbstractString)
    open(kernelfile) do f
        return parse_kernel(f)
    end
end

function parse_kernel(stream::IO)

    sections = Vector{Section}()
    #current = Section("__julia__", Vector{AbstractString}())
    current = Section("__julia__", AbstractString[])
    push!(sections, current)
    
    for line in eachline(stream)
        s = strip(line)

        if length(s) > 2 && s[1] == '[' && s[end] == ']'
            current = Section(s[2:end-1], Vector{AbstractString}())
            push!(sections, current)

        else
            push!(current.body, line)

        end
    end
    sections
end

end
