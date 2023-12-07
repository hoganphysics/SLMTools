module SubSubmodule
using ..Submodule: subtestA, subtestB
export subsubtestA, subsubtestB

function subsubtestA()
    println("subsubtestA")
end

function subsubtestB()
    println("subsubtestB but also")
    subtestA()
end

end # module SubSubmodule
