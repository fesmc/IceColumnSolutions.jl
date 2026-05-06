module IceColumnSolutions

include("params.jl")
include("solution.jl")
include("stationary.jl")
include("transient.jl")
include("benchmarks.jl")

export IceColumnPar, IceColumn
export solve_stationary, solve
export uniform, stationary_init, eigenvalues
export benchmark
export to_celsius, to_celsius!

end
