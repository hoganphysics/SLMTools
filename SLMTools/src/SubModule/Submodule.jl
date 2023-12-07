module Submodule

export subtestA, subtestB

include("SubSubmodule.jl") #wont work here
# module SubSubmodule
# using ..Submodule: subtestA, subtestB
# export subsubtestA, subsubtestB

# function subsubtestA()
#     println("subsubtestA")
# end

# function subsubtestB()
#     println("subsubtestB but also")
#     subtestA()
# end

# end

using .SubSubmodule
export subsubtestA
export subsubtestB


function subtestA()
    println("subtestA")
end

function subtestB()
    println("subtestB")
end


end # module Submodule
